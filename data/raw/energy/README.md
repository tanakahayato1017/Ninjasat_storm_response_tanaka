# NinjaSat specific orbital energy (E_J2) — pending release

`specific_energy.csv` holds the $J_2$-corrected specific orbital energy
$E_{J2}$ computed at each GNSS epoch (native ~10-min cadence, full 2024)
directly from the GNSS-derived position and velocity:

$$E_{J2} = \tfrac{1}{2}v_i^2 - \frac{\mu}{r}\left[1 - \frac{J_2}{2}\left(\frac{R_E}{r}\right)^2(3\sin^2\phi - 1)\right]$$

This quantity is treated the same as raw GNSS telemetry for release
purposes: although it is a derived scalar rather than a raw position fix,
it is close enough to a direct restatement of orbital radius (hence
altitude) at each epoch that its release is pending the NinjaSat team's
confirmation, on the same terms as [`../GNSS_data/`](../GNSS_data/).

`figures/fig1_lomb_scargle.py` reads from this path and will run once the
file is added here.
