# Redacted energy-rate data

`dEdt_storm_3h_redacted.csv` is a public-safe copy of the internal
`dEdt_storm_3h.csv` product. It retains only `datetime`, `event`, and
`dEdt_Jkg_per_day`, which are required by the Joule-heating and density
validation scripts.

The source file's absolute `E_J2_mean` specific orbital energy was removed
because it can be inverted to an absolute orbital radius and is pending
NinjaSat-team release. Other unused columns were also omitted. Neither
consumer used any removed column in its original numerical calculation.
