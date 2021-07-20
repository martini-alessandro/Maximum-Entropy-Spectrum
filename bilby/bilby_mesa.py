#!/usr/bin/env python
"""
A script to run bilby with a PSD calculated by MESA on the residual

Usage
-----
To run the script, first ensure both bilby and pymultinest are installed

    $ conda install -c conda-forge bilby pymultinest

Then to run using the MESA likelihood:

    $ python bilby_mesa.py --model mesa --plot-psd

Other options can be found from

    $ python bilby_mesa.py --help

"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
import bilby
from bilby.gw.likelihood import Likelihood
import memspectrum
from scipy.signal.windows import tukey


parser = argparse.ArgumentParser(__doc__)
parser.add_argument("-s", "--seed", default=42, type=int)
parser.add_argument("-m", "--model", default="mesa")
parser.add_argument("--nlive", default=1000)
parser.add_argument("--plot-psd", action="store_true")
parser.add_argument("--mesa-m", default=None, type=int)
parser.add_argument("--mesa-optimisation-method", default="FPE")
parser.add_argument("--mesa-method", default="Fast")
args, _ = parser.parse_known_args()


class MESAGravitationalWaveTransient(Likelihood):
    def __init__(
        self,
        interferometers,
        waveform_generator,
        mesa_m,
        mesa_method,
        mesa_optimisation_method,
    ):
        super(MESAGravitationalWaveTransient, self).__init__(dict())

        self.interferometers = interferometers
        self.waveform_generator = waveform_generator
        self._nll = np.nan
        self.mesa = memspectrum.MESA()
        self.mesa_m = mesa_m
        self.mesa_method = mesa_method
        self.mesa_optimisation_method = mesa_optimisation_method

    def update_ifo_psd(self, interferometer, htilde=None):
        """ Use MESA to update the PSD for the given interferometer

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
        alpha = 2 * roll_off / duration
        window = tukey(N, alpha=alpha)
        delta*= window

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

    def log_likelihood(self):
        # Calculate the waveform_polarizations dictionary
        waveform_polarizations = self.waveform_generator.frequency_domain_strain(
            self.parameters.copy()
        )

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
        termB = -0.5 * np.sum(np.log(0.5 * np.pi * duration * psd))

        log_l = termA + termB
        return log_l.real

    def noise_log_likelihood(self):
        log_l = 0
        for interferometer in self.interferometers:
            log_l += self.noise_log_likelihood_interferometer(interferometer)
        return log_l.real

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
        termA = (
            -2.0
            * np.sum(abs(data) ** 2 / psd)
            / self.waveform_generator.duration
        )

        # Calculate the second term in Eq (8) of Veitch (2015)
        termB = -0.5 * np.sum(np.log(0.5 * np.pi * duration * psd))

        log_l = termA + termB
        return log_l.real


# Set up of the simulation
duration = 4.0
sampling_frequency = 1024.0
maximum_frequency = sampling_frequency / 2
minimum_frequency = 20
minimum_frequency_simulation = 10
outdir = "outdir"
np.random.seed(args.seed)

# Hold these parameters fixed for the simulation
fixed_parameters = dict(
    a_1=0.0,
    a_2=0.0,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    luminosity_distance=2000.0,
    theta_jn=0.4,
    psi=2.659,
    phase=1.3,
    geocent_time=1126259642.413,
    ra=1.375,
    dec=-1.2108,
)

priors = bilby.gw.prior.BBHPriorDict()
priors["geocent_time"] = bilby.core.prior.Uniform(
    minimum=fixed_parameters["geocent_time"] - 0.1,
    maximum=fixed_parameters["geocent_time"] + 0.1,
    name="geocent_time",
    latex_label="$t_c$",
    unit="$s$",
)
priors["chirp_mass"] = bilby.core.prior.Uniform(
    minimum=35, maximum=40, name="chirp_mass"
)
for key, val in fixed_parameters.items():
    priors[key] = bilby.core.prior.DeltaFunction(val)

injection_parameters = priors.sample()

waveform_arguments = dict(
    waveform_approximant="TaylorF2",
    reference_frequency=minimum_frequency,
    minimum_frequency=minimum_frequency,
    catch_waveform_errors=True,
)
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)

ifos = bilby.gw.detector.InterferometerList(["H1"])
for ifo in ifos:
    ifo.minimum_frequency = minimum_frequency_simulation
    ifo.maximum_frequency = maximum_frequency
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=fixed_parameters["geocent_time"] + 2 - duration,
)
ifos.inject_signal(
    waveform_generator=waveform_generator, parameters=injection_parameters
)

# Store the simulation PSDs for plotting later on
true_psds = [ifo.power_spectral_density_array for ifo in ifos]

# Set of arguments for the sampler
skwargs = dict(sampler="pymultinest", nlive=int(args.nlive), dlogz=0.5)

# A label to store the results by
label = f"{args.model}_seed{args.seed}"
if args.model == "mesa":
    label += f"_m{args.mesa_m}_{args.mesa_method}_{args.mesa_optimisation_method}"

if args.model == "standard":
    standard_likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers=ifos, waveform_generator=waveform_generator
    )

    result = bilby.run_sampler(
        likelihood=standard_likelihood,
        priors=priors,
        injection_parameters=injection_parameters,
        outdir=outdir,
        label=label,
        use_ratio=True,
        **skwargs,
    )
else:
    mesa_likelihood = MESAGravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=waveform_generator,
        mesa_m=args.mesa_m,
        mesa_method=args.mesa_method,
        mesa_optimisation_method=args.mesa_optimisation_method,
    )

    result = bilby.run_sampler(
        likelihood=mesa_likelihood,
        priors=priors,
        injection_parameters=injection_parameters,
        outdir=outdir,
        label=label,
        use_ratio=True,
        **skwargs,
    )

    if args.plot_psd:
        print("Generating PSD plot")
        for i, ifo in enumerate(ifos):
            fig, ax = plt.subplots()
            ax.loglog(
                ifo.frequency_array, np.abs(ifo.frequency_domain_strain)**2,
                color="C0", label="Data"
            )
            for _ in range(500):
                sample = dict(result.posterior.sample().iloc[0])
                waveform_polarizations = waveform_generator.frequency_domain_strain(
                    sample
                )
                htilde = ifo.get_detector_response(waveform_polarizations, sample)
                mesa_likelihood.parameters.update(sample)
                ifo = mesa_likelihood.update_ifo_psd(ifo, htilde)
                ax.loglog(
                    ifo.frequency_array,
                    ifo.power_spectral_density_array,
                    color="tab:red",
                    zorder=100,
                    alpha=0.8,
                )
            ax.loglog(ifo.frequency_array, true_psds[i], color="C1", alpha=0.8, label="True PSD")
            ax.set_xlim(minimum_frequency - 10, maximum_frequency + 100)
            ax.axvspan(minimum_frequency, maximum_frequency, color="k", alpha=0.1)
            ax.set_xlabel("Frequency [Hz]")
            ax.set_ylabel("PSD [1/Hz]")
            ax.legend()
            fig.savefig(f"outdir/{label}_{ifo.name}_psd.png", dpi=500)
            fig.clf()

result.plot_corner()
