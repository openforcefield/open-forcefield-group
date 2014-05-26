# Density data for neat liquids and binary mixtures

## Summary

This repository contains temperature-dependent density data at ambient pressure for neat liquids and binary mixtures.

## Protocol

All measurements were conducted with a [Mettler-Toledo DM40 LiquiPhysics Excellence density meter](http://us.mt.com/us/en/home/products/Laboratory_Analytics_Browse/Density_Family_Browse_main/DE_Benchtop/DM40.html).
This meter has an accuracy of 0.0001 g/cm3 and a range of 0-3 g/cm3.
Sample volume requirements are approximately 1 mL.
The instrument can read in a temperature range of 0-91 C, requiring approximately 30 s per measurement.

A [Mettler-Toledo SC30 sample delivery and cleaning unit](http://us.mt.com/us/en/home/products/Laboratory_Analytics_Browse/Density_Family_Browse_main/DE_autom_sys/SC30.html) with 30-sample capacity was used to automate the collection of data at different mole fractions of composition for binary mixtures.

All samples will be prepared using a [Mettler-Toledo Quantos automated liquid dosing system](http://us.mt.com/us/en/home/products/Laboratory_Weighing_Solutions/aut_do_sys_qua/liquid_dosing.html) to dispense neat liquids or binary mixtures under argon.
Binary mixtures will be dispensed at 30 different mole fractions, equally spaced from x = 0 to 1, though actual mole fractions will differ slightly due to dispensing inaccuracies.
The balance is a [Mettler-Toledo QD206DR](http://us.mt.com/us/en/home/products/Laboratory_Weighing_Solutions/aut_do_sys_qua/liquid_dosing.html), with readability (81 g x 0.005 mg / 220 g x 0.01 mg).

Define:
* `x`, `dx` - mole fraction and its uncertainty
* `a`, `b` - moles of liquids dispensed
* `da`, `db` - uncertainty in moles of dispensed liquids

Note that the moles of liquid dispensed is equal to the mass of liquid dispensed divided by its molecular weight.

We have
```
x = a / (a + b)
dx = x * sqrt[ (da/a)^2 + (db/b)^2 ]
```
For a 1 mL ~ 1 g sample, we estimate the error in x (for x ~ 0.5) should be roughly
```
dx / x = sqrt[ (da/a)^2 + (db/b)^2 ] ~ sqrt[ (0.01 / 500)^2 + (0.01 / 500)^2 ] ~ 0.00002
```
or 0.002%.

## Useful references

* A database of molecules and liquids with their CRC data and GAFF/OPLS parameter files: http://virtualchemistry.org

