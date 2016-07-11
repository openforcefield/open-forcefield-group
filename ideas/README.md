# Ideas for forcefield parameterization techniques

## Manifest
* `bayesian-atomtype-sampling` - Illustration of Bayesian sampling over atom SMARTS-based atomtypes
* `bayesian-gbsa-parameterization` - Illustration of Bayesian GBSA parameterization using single-point small molecule hydration free energies

## Glossary
* Smarty - SMARTS-based chemical perception of atom types, used as a test of the SMARTS-based chemical perception for force fields.  Involves the Bayesian sampling of SMARTS (or SMIRKS) strings to represent atom types rather than unique string assigned by intuition or experience.
* SMIRFF - SMIRKS-based force field. Force fields are assigned based on SMIRKS strings determined by the chemical environments, rather than expert-devised atom typing. The key is going directly from the chemical perception to the force field, instead of through the intermediate step of atom type.
* BaFFLE - Bayesian sampling of chemical perception or force field parameters to produce a force field.
