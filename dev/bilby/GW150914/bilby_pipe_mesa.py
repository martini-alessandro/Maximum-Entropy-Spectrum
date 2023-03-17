#!/usr/bin/env python
import sys
import signal

from bilby_pipe.data_analysis import DataAnalysisInput, sighandler
from bilby_pipe.parser import create_parser
from bilby_pipe.main import parse_args
from bilby_pipe.utils import (
    CHECKPOINT_EXIT_CODE, log_version_information, logger
)

import bilby
import numpy as np
import memspectrum
from scipy.signal.windows import tukey


class MESAGravitationalWaveTransient(bilby.gw.likelihood.GravitationalWaveTransient):
    def __init__(
        self,
        interferometers,
        waveform_generator,
        mesa_m,
        mesa_method,
        mesa_optimisation_method,
        reference_frame="sky",
    ):
        super(MESAGravitationalWaveTransient, self).__init__(
            interferometers=interferometers,
            waveform_generator=waveform_generator,
            reference_frame=reference_frame,
        )

        self.interferometers = interferometers
        self.waveform_generator = waveform_generator
        self._nll = np.nan
        self.mesa = memspectrum.MESA()
        self.mesa_m = mesa_m
        self.mesa_method = mesa_method
        self.mesa_optimisation_method = mesa_optimisation_method
        self._cached_noise_log_likelihood_dict = {}

    def update_ifo_psd(self, interferometer, htilde=None):
        """Use MESA to update the PSD for the given interferometer

        Parameters
        ----------
        interferometer: bilby.gw.interferometer.Interferometer
            The bilby interferometer object which stores its own PSD.
        htilde: np.array
            The complex frequency-domain detector response to the signal.
            If this is None, no signal is subtracted from the time-domain
            data before calculated the PSD.

        """
        if htilde is not None:
            h_timedomain = bilby.core.utils.infft(
                htilde, self.waveform_generator.sampling_frequency
            )
        else:
            h_timedomain = 0

        delta = interferometer.time_domain_strain - h_timedomain

        N = len(delta)
        roll_off = 0.2  # Rise time of window in seconds
        alpha = 2 * roll_off / interferometer.duration
        window = tukey(N, alpha=alpha)
        delta *= window

        _ = self.mesa.solve(
            delta,
            m=self.mesa_m,
            method=self.mesa_method,
            optimisation_method=self.mesa_optimisation_method,
        )
        _, psd = self.mesa.spectrum(
            dt=1 / interferometer.sampling_frequency,
            onesided=True,
        )

        cached_psd = interferometer.power_spectral_density._cache["psd_array"]
        psd = np.concatenate([psd, [cached_psd[-1]]])
        psd[~interferometer.frequency_mask] = np.inf
        interferometer.power_spectral_density._cache["psd_array"] = psd
        return interferometer

    def log_likelihood_ratio(self):
        return self.log_likelihood() - self.noise_log_likelihood()

    def log_likelihood(self):

        # Update m
        self.mesa_m = int(self.parameters["mesa_m"])

        # Effective prior cut on m=0
        if self.mesa_m == 0:
            return -np.nan_to_num(np.inf)

        # Calculate the waveform_polarizations dictionary
        waveform_polarizations = self.waveform_generator.frequency_domain_strain(
            self.parameters.copy()
        )

        self.parameters.update(self.get_sky_frame_parameters())

        # If the waveform_polarizations dict is None, calculation failed
        # so return -inf likelihood (failure usualy means we outstepped the
        # bounds of the waveform model.
        if waveform_polarizations is None:
            return np.nan_to_num(-np.inf)

        # Sum up the log_likelihood per detector
        log_l = 0
        for interferometer in self.interferometers:
            log_l += self.log_likelihood_interferometer(
                waveform_polarizations, interferometer
            )

        return log_l.real

    def log_likelihood_interferometer(self, waveform_polarizations, interferometer):
        """

        Parameters
        ==========
        waveform_polarizations: dict
            Dictionary containing the desired waveform polarization modes and the related strain
        interferometer: bilby.gw.detector.Interferometer
            The Interferometer object we want to have the log-likelihood for

        Returns
        =======
        float: The real part of the log-likelihood for this interferometer

        """
        # Get a dictionary of the + and x response
        htilde = interferometer.get_detector_response(
            waveform_polarizations, self.parameters
        )

        # Update the interferometer PSD using MESA
        interferometer = self.update_ifo_psd(interferometer, htilde)

        # Extract the PSD and frequency_domain_strain data
        psd = interferometer.power_spectral_density_array
        data = interferometer.frequency_domain_strain

        # Apply the frequency mask (cutting data outdir minimum_frequency-maximum_frequency)
        htilde = htilde[interferometer.frequency_mask]
        psd = psd[interferometer.frequency_mask]
        data = data[interferometer.frequency_mask]

        # Calculate the first term in Eq (10) of Veitch (2015)
        termA = (
            -2.0
            * np.vdot(data - htilde, (data - htilde) / psd)
            / self.waveform_generator.duration
        )

        # Calculate the second term in Eq (10) of Veitch (2015)
        termB = -0.5 * np.sum(np.log(0.5 * np.pi * interferometer.duration * psd))

        log_l = termA + termB
        return log_l.real

    def noise_log_likelihood(self):
        self.mesa_m = int(self.parameters["mesa_m"])

        # Effective prior cut on m=0
        if self.mesa_m == 0:
            return -np.nan_to_num(np.inf)

        if self.mesa_m in self._cached_noise_log_likelihood_dict:
            return self._cached_noise_log_likelihood_dict[self.mesa_m]
        log_l = 0
        for interferometer in self.interferometers:
            log_l += self.noise_log_likelihood_interferometer(interferometer)
        noise_log_likelihood = log_l.real
        self._cached_noise_log_likelihood_dict[self.mesa_m] = noise_log_likelihood
        return noise_log_likelihood

    def noise_log_likelihood_interferometer(self, interferometer):
        # Update the interferometer PSD using MESA
        interferometer = self.update_ifo_psd(interferometer)

        # Extract the PSD and frequency_domain_strain data
        psd = interferometer.power_spectral_density_array
        data = interferometer.frequency_domain_strain

        # Apply the frequency mask (cutting data outdir minimum_frequency-maximum_frequency)
        psd = psd[interferometer.frequency_mask]
        data = data[interferometer.frequency_mask]

        # Calculate the first term in Eq (8) of Veitch (2015)
        termA = -2.0 * np.sum(abs(data) ** 2 / psd) / self.waveform_generator.duration

        # Calculate the second term in Eq (8) of Veitch (2015)
        termB = -0.5 * np.sum(np.log(0.5 * np.pi * interferometer.duration * psd))

        log_l = termA + termB
        return log_l.real


def main():
    """ Data analysis main logic """
    args, unknown_args = parse_args(sys.argv[1:], create_parser(top_level=False))
    log_version_information()
    analysis = DataAnalysisInput(args, unknown_args)

    # Get the standard likelihood and priors
    standard_likelihood, run_priors = analysis.get_likelihood_and_priors()
    likelihood = MESAGravitationalWaveTransient(
        interferometers=standard_likelihood.interferometers,
        waveform_generator=standard_likelihood.waveform_generator,
        mesa_m=10,
        mesa_method="Fast",
        mesa_optimisation_method="Fixed",
        reference_frame=standard_likelihood.reference_frame
    )

    #prop_priors = bilby.gw.prior.BBHPriorDict(run_priors)
    #prop_priors.pop('mesa_m')
    #plist = bilby.bilby_mcmc.proposals.get_proposal_cycle("gwA", prop_priors).proposal_list
    #plist.append(
        #bilby.bilby_mcmc.proposals.PriorProposal(
            #run_priors, subset=["mesa_m"], weight=10
        #)
    #)
    #analysis.sampler_kwargs["proposal_cycle"] = bilby.bilby_mcmc.proposals.ProposalCycle(plist)

    if analysis.scheduler.lower() == "condor":
        signal.signal(signal.SIGALRM, handler=sighandler)
        signal.alarm(analysis.periodic_restart_time)

    bilby.run_sampler(
        likelihood=likelihood,
        priors=run_priors,
        sampler=analysis.sampler,
        label=analysis.label,
        outdir=analysis.result_directory,
        conversion_function=analysis.parameter_generation,
        injection_parameters=analysis.meta_data["injection_parameters"],
        meta_data=analysis.meta_data,
        result_class=analysis.result_class,
        exit_code=CHECKPOINT_EXIT_CODE,
        save=analysis.result_format,
        **analysis.sampler_kwargs,
    )

    sys.exit(0)


if __name__ == "__main__":
    main()
