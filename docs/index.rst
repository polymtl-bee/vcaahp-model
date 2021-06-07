Model variable capacity air-to-air heat pumps in TRNSYS
===========================================================

The TRNSYS Types 3254 and 3223 presented here aim to model single zone,
mini-split air-to-air heat pumps equipped with a variable speed comp­ressor.
The Type 3254 relies on :ref:`performance maps <performance file>` to
model heat pump performance, and the Type 3223 is dedicated to control aspects
(see figure below). As of this version, the model is exclusively dedicated
to short time step simulations—ideally below 5 minutes.

.. figure:: pictures/model-structure.pdf

   Structure of the VCAAHP model.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   type3223/type3223
   type3254/type3254


Features
--------

- Heating and cooling mode operation
- Compressor frequency modulation with a dedicated controller
- Defrost cycles in heating mode
- Sensible and latent cooling


Support
-------

If you encounter bugs or are having problems, please open an issue in the
`issue tracker <https://github.com/polymtl-bee/vcaahp-model/issues>`_
and submit a minimal working example in the form of a deck file
to highlight what is not working.

You can ask questions about the model in the Q&A Discussions.
