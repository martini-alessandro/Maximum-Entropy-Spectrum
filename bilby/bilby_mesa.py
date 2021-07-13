#!/usr/bin/env python
import sys

import numpy as np
import matplotlib.pyplot as plt
import bilby
from bilby.gw.likelihood import Likelihood
import memspectrum


class MESAGravitationalWaveTransient(Likelihood):
    def __init__(self, interferometers, waveform_generator):
        super(MESAGravitationalWaveTransient, self).__init__(dict())

        self.interferometers = interferometers
        self.waveform_generator = waveform_generator
        self._nll = np.nan
        self.m = memspectrum.MESA()

    def log_likelihood(self):
        log_l = 0
        waveform_polarizations =\
            self.waveform_generator.frequency_domain_strain(
                self.parameters.copy())
        if waveform_polarizations is None:
            return np.nan_to_num(-np.inf)
        for interferometer in self.interferometers:
            log_l += self.log_likelihood_interferometer(
                waveform_polarizations, interferometer)
        return log_l.real

    def update_ifo_psd(self, interferometer, signal_ifo):
        signal_projected = bilby.core.utils.infft(
            signal_ifo,
            self.waveform_generator.sampling_frequency
        )

        delta = interferometer.time_domain_strain - signal_projected
        _ = self.m.solve(delta, m=100)
        _, psd = self.m.spectrum(
            dt=1 / interferometer.sampling_frequency,
            onesided=True,
        )

        psd = np.concatenate([psd, [np.inf]])
        interferometer.power_spectral_density._cache["psd_array"] = psd
        return interferometer

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

        signal_ifo = interferometer.get_detector_response(
            waveform_polarizations, self.parameters)

        interferometer = self.update_ifo_psd(interferometer, signal_ifo)

        log_l = - 2. / self.waveform_generator.duration * np.vdot(
            interferometer.frequency_domain_strain - signal_ifo,
            (interferometer.frequency_domain_strain - signal_ifo) /
            interferometer.power_spectral_density_array)

        log_l += 2. / self.waveform_generator.duration * np.sum(
            abs(interferometer.frequency_domain_strain) ** 2 /
            interferometer.power_spectral_density_array)

        return log_l.real


duration = 4.0
sampling_frequency = 512.0
outdir = "outdir"
np.random.seed(150914)

injection_parameters = dict(
    chirp_mass=36.0,
    mass_ratio=0.9,
    a_1=0.0,
    a_2=0.0,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    luminosity_distance=4000.0,
    theta_jn=0.4,
    psi=2.659,
    phase=1.3,
    geocent_time=1126259642.413,
    ra=1.375,
    dec=-1.2108,
)

waveform_arguments = dict(
    waveform_approximant="TaylorF2",
    reference_frequency=20.0,
    minimum_frequency=20.0,
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
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 3,
)
ifos.inject_signal(
    waveform_generator=waveform_generator, parameters=injection_parameters
)

priors = bilby.gw.prior.BBHPriorDict()
priors["geocent_time"] = bilby.core.prior.Uniform(
    minimum=injection_parameters["geocent_time"] - 0.1,
    maximum=injection_parameters["geocent_time"] + 0.1,
    name="geocent_time",
    latex_label="$t_c$",
    unit="$s$",
)
priors["chirp_mass"] = bilby.core.prior.Uniform(
    minimum=35, maximum=37, name="chirp_mass"
)

for key in bilby.gw.source.spin:
    if key in injection_parameters:
        priors[key] = injection_parameters[key]
for key in [
    "psi",
    "ra",
    "dec",
    "geocent_time",
    "phase",
    "theta_jn",
    "luminosity_distance",
]:
    priors[key] = injection_parameters[key]


true_psds = [ifo.power_spectral_density_array for ifo in ifos]

skwargs = dict(sampler="pymultinest", nlive=1000, dlogz=0.5)

if "standard" in sys.argv:
    standard_likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers=ifos, waveform_generator=waveform_generator
    )

    result = bilby.run_sampler(
        likelihood=standard_likelihood,
        priors=priors,
        injection_parameters=injection_parameters,
        outdir=outdir,
        label="standard",
        use_ratio=True,
        **skwargs
    )
else:
    mesa_likelihood = MESAGravitationalWaveTransient(
        interferometers=ifos, waveform_generator=waveform_generator
    )

    result = bilby.run_sampler(
        likelihood=mesa_likelihood,
        priors=priors,
        injection_parameters=injection_parameters,
        outdir=outdir,
        label="mesa",
        use_ratio=False,
        **skwargs
    )

    for i, ifo in enumerate(ifos):
        fig, ax = plt.subplots()
        print("Generating PSD plot")
        for _ in range(1000):
            sample = dict(result.posterior.sample().iloc[0])
            mesa_likelihood.parameters.update(sample)
            waveform_polarizations =\
                waveform_generator.frequency_domain_strain(sample)
            signal_ifo = ifo.get_detector_response(
                waveform_polarizations, sample)
            ifo = mesa_likelihood.update_ifo_psd(ifo, signal_ifo)
            ax.loglog(ifo.frequency_array, ifo.power_spectral_density_array, color='k', lw=0.5)
        ax.loglog(ifo.frequency_array, true_psds[i], ls='--')
        ax.set_xlim(15, 256)
        ax.set_xlabel("Frequency [Hz]")
        fig.savefig(f"{ifo.name}_psd.png", dpi=500)

result.plot_corner()
