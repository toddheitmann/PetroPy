.. image:: https://toddheitmann.github.io/PetroPy/_images/petropy_logo.png

PetroPy
=======

A python petrophysics package allowing scientific python computing
of conventional and unconventional formation evaluation. Reads las
files using `lasio <https://github.com/kinverarity1/lasio>`__. Includes
a petrophysical workflow and a log viewer based on XML templates.

.. image:: https://toddheitmann.github.io/PetroPy/_images/university_6-18W_no1.png

Requirements
------------

-  `cchardet <https://github.com/PyYoshi/cChardet>`__
-  `lasio <https://github.com/kinverarity1/lasio>`__
-  `numpy <http://www.numpy.org>`__
-  `scipy <https://www.scipy.org>`__
-  `pandas <http://pandas.pydata.org>`__
-  `matplotlib <http://matplotlib.org>`__
-  `scikit-learn <http://scikit-learn.org/stable/>`__

Installation
------------

Install PetroPy through pip via the command line

.. code-block:: bash

  pip install petropy

To read in an las file, pass the file reference:

.. code-block:: python

  import petropy as ptr
  file_path = r'path/to/well.las'
  log = ptr.Log(file_path)

Documentation
-------------

View the `online documentation`_ for classes and methods.

.. _online documentation: https://toddheitmann.github.io/PetroPy/

Las File Processing
-------------------

To understanding using petropy in a petrophysical workflow for las file
processing, see the `example page`_.

.. _example page: https://toddheitmann.github.io/PetroPy/auto_examples/

Petrophysical Model Quick Look
------------------------------

.. code-block:: python

  >>> # import petropy and print raw curves
  >>> import petropy as ptr
  >>> log = ptr.log_data('WFMP')
  >>> print(log.curves)

.. code-block:: bash

  Mnemonic   Unit  Value         Description
  --------   ----  -----         -----------
  DEPT       F     00 000 00 00  1  Depth Curve
  CALI       INCH  99 075 22 05  2  CALIPER
  DPHI       DECP  99 075 22 05  3  DENSITY POROSITY -LIME-
  GR         GAPI  99 075 22 05  4  GAMMA RAY
  NPHI       DECP  99 075 22 05  5  NEUTRON POROSITY -LIME-
  PE         B/E   99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB       G/C3  99 075 22 05  7  BULK DENSITY
  PHIX       DECP  99 075 22 05  8  CROSSPLOT POROSITY
  C13        INCH  99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  C24        INCH  99 075 22 05  10  CALIPER PADS 2 - 4    -FACT-
  DT         US/F  99 075 22 05  11  SONIC TRANSIT TIME
  SPHI       DECP  99 075 22 05  12  SONIC POROSITY  -LIME-
  GR3              99 075 22 05  13  GAMMA RAY
  ILD        OHMM  99 075 22 05  14  IL, DEEP RESISTIVITY
  ILM        OHMM  99 075 22 05  15  IL, MEDIUM RESISTIVITY
  SGRD       OHMM  99 075 22 05  16  SHORT GUARD RESISTIVITY
  SP         MV    99 075 22 05  17  SPONTANEOUS POTENTIAL
  CAL_N      INCH  99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  GR_N       GAPI  99 075 22 05  4  GAMMA RAY
  RESMED_N   OHMM  99 075 22 05  15  IL, MEDIUM RESISTIVITY
  RESDEEP_N  OHMM  99 075 22 05  14  IL, DEEP RESISTIVITY
  NPHI_N     DECP  99 075 22 05  5  NEUTRON POROSITY -LIME-
  DPHI_N     DECP  99 075 22 05  3  DENSITY POROSITY -LIME-
  SPHI_N     DECP  99 075 22 05  12  SONIC POROSITY  -LIME-
  PE_N       B/E   99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB_N     G/C3  99 075 22 05  7  BULK DENSITY
  DTC_N      US/F  99 075 22 05  11  SONIC TRANSIT TIME
  SP_N       MV    99 075 22 05  17  SPONTANEOUS POTENTIAL

.. code-block:: python

  >>> # read tops into Log object and print
  >>> log.tops_from_csv()
  >>> print(log.tops)

.. code-block:: bash

  {'WFMPA': 6993.5, 'WFMPB': 7294.0, 'WFMPC': 7690.5, 'WFMPD': 8028.0}

.. code-block:: python

  >>> # load default parameters and print values
  >>> log.fluid_properties_parameters_from_csv()
  >>> print(log.fluid_properties_parameters.keys())

.. code-block:: bash

  dict_keys(['default', 'WFMP'])

.. code-block:: python

  >>> # specificy formation intervals
  >>> f = ['WFMPA', 'WFMPB', 'WFMPC']
  >>> # calculate fluid properties for defined formations
  >>> log.formation_fluid_properties(f, parameter = 'WFMP')
  >>> # print curves for description of calculated curves
  >>> print(log.curves)

.. code-block:: bash

  Mnemonic    Unit  Value         Description
  --------    ----  -----         -----------
  DEPT        F     00 000 00 00  1  Depth Curve
  CALI        INCH  99 075 22 05  2  CALIPER
  DPHI        DECP  99 075 22 05  3  DENSITY POROSITY -LIME-
  GR          GAPI  99 075 22 05  4  GAMMA RAY
  NPHI        DECP  99 075 22 05  5  NEUTRON POROSITY -LIME-
  PE          B/E   99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB        G/C3  99 075 22 05  7  BULK DENSITY
  PHIX        DECP  99 075 22 05  8  CROSSPLOT POROSITY
  C13         INCH  99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  C24         INCH  99 075 22 05  10  CALIPER PADS 2 - 4    -FACT-
  DT          US/F  99 075 22 05  11  SONIC TRANSIT TIME
  SPHI        DECP  99 075 22 05  12  SONIC POROSITY  -LIME-
  GR3               99 075 22 05  13  GAMMA RAY
  ILD         OHMM  99 075 22 05  14  IL, DEEP RESISTIVITY
  ILM         OHMM  99 075 22 05  15  IL, MEDIUM RESISTIVITY
  SGRD        OHMM  99 075 22 05  16  SHORT GUARD RESISTIVITY
  SP          MV    99 075 22 05  17  SPONTANEOUS POTENTIAL
  CAL_N       INCH  99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  GR_N        GAPI  99 075 22 05  4  GAMMA RAY
  RESMED_N    OHMM  99 075 22 05  15  IL, MEDIUM RESISTIVITY
  RESDEEP_N   OHMM  99 075 22 05  14  IL, DEEP RESISTIVITY
  NPHI_N      DECP  99 075 22 05  5  NEUTRON POROSITY -LIME-
  DPHI_N      DECP  99 075 22 05  3  DENSITY POROSITY -LIME-
  SPHI_N      DECP  99 075 22 05  12  SONIC POROSITY  -LIME-
  PE_N        B/E   99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB_N      G/C3  99 075 22 05  7  BULK DENSITY
  DTC_N       US/F  99 075 22 05  11  SONIC TRANSIT TIME
  SP_N        MV    99 075 22 05  17  SPONTANEOUS POTENTIAL
  PORE_PRESS  psi                 Calculated Pore Pressure
  RES_TEMP    F                   Calculated Reservoir Temperature
  NES         psi                 Calculated Net Effective Stress
  RW          ohmm                Calculated Resistivity Water
  RMF         ohmm                Calculated Resistivity Mud Filtrate
  RHO_HC      g/cc                Calculated Density of Hydrocarbon
  RHO_W       g/cc                Calculated Density of Water
  RHO_MF      g/cc                Calculated Density of Mud Filtrate
  NPHI_HC     v/v                 Calculated Neutron Log Response of Hydrocarbon
  NPHI_W      v/v                 Calculated Neutron Log Response of Water
  NPHI_MF     v/v                 Calculated Neutron Log Response of Mud Filtrate
  MU_HC       cP                  Calculated Viscosity of Hydrocarbon
  BO                              Calculated Oil Formation Volume Factor
  BP          psi                 Calcualted Bubble Point

.. code-block:: python

  >>> # load default multimineral parameters
  >>> log.multimineral_parameters_from_csv()
  >>> # print available default formation parameters
  >>> print(log.multimineral_parameters.keys())

.. code-block:: bash

  dict_keys(['default', 'WFMP'])

.. code-block:: python

  >>> # calculate mulitmineral model over defined formations
  >>> # with parameter 'WFMP'
  >>> log.formation_multimineral_model(f, parameter = 'WFMP')
  >>> log.write('processed_log.las')
  >>> # print curves for description of calculated curves
  >>> print(log.curves)

.. code-block:: bash

  Mnemonic    Unit   Value         Description
  --------    ----   -----         -----------
  DEPT        F      00 000 00 00  1  Depth Curve
  CALI        INCH   99 075 22 05  2  CALIPER
  DPHI        DECP   99 075 22 05  3  DENSITY POROSITY -LIME-
  GR          GAPI   99 075 22 05  4  GAMMA RAY
  NPHI        DECP   99 075 22 05  5  NEUTRON POROSITY -LIME-
  PE          B/E    99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB        G/C3   99 075 22 05  7  BULK DENSITY
  PHIX        DECP   99 075 22 05  8  CROSSPLOT POROSITY
  C13         INCH   99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  C24         INCH   99 075 22 05  10  CALIPER PADS 2 - 4    -FACT-
  DT          US/F   99 075 22 05  11  SONIC TRANSIT TIME
  SPHI        DECP   99 075 22 05  12  SONIC POROSITY  -LIME-
  GR3                99 075 22 05  13  GAMMA RAY
  ILD         OHMM   99 075 22 05  14  IL, DEEP RESISTIVITY
  ILM         OHMM   99 075 22 05  15  IL, MEDIUM RESISTIVITY
  SGRD        OHMM   99 075 22 05  16  SHORT GUARD RESISTIVITY
  SP          MV     99 075 22 05  17  SPONTANEOUS POTENTIAL
  CAL_N       INCH   99 075 22 05  9  CALIPER PADS 1 - 3    -FACT-
  GR_N        GAPI   99 075 22 05  4  GAMMA RAY
  RESMED_N    OHMM   99 075 22 05  15  IL, MEDIUM RESISTIVITY
  RESDEEP_N   OHMM   99 075 22 05  14  IL, DEEP RESISTIVITY
  NPHI_N      DECP   99 075 22 05  5  NEUTRON POROSITY -LIME-
  DPHI_N      DECP   99 075 22 05  3  DENSITY POROSITY -LIME-
  SPHI_N      DECP   99 075 22 05  12  SONIC POROSITY  -LIME-
  PE_N        B/E    99 075 22 05  6  PHOTO-ELECTRIC FACTOR
  RHOB_N      G/C3   99 075 22 05  7  BULK DENSITY
  DTC_N       US/F   99 075 22 05  11  SONIC TRANSIT TIME
  SP_N        MV     99 075 22 05  17  SPONTANEOUS POTENTIAL
  PORE_PRESS  psi                  Calculated Pore Pressure
  RES_TEMP    F                    Calculated Reservoir Temperature
  NES         psi                  Calculated Net Effective Stress
  RW          ohmm                 Calculated Resistivity Water
  RMF         ohmm                 Calculated Resistivity Mud Filtrate
  RHO_HC      g/cc                 Calculated Density of Hydrocarbon
  RHO_W       g/cc                 Calculated Density of Water
  RHO_MF      g/cc                 Calculated Density of Mud Filtrate
  NPHI_HC     v/v                  Calculated Neutron Log Response of Hydrocarbon
  NPHI_W      v/v                  Calculated Neutron Log Response of Water
  NPHI_MF     v/v                  Calculated Neutron Log Response of Mud Filtrate
  MU_HC       cP                   Calculated Viscosity of Hydrocarbon
  BO                               Calculated Oil Formation Volume Factor
  BP          psi                  Calcualted Bubble Point
  PHIE        v/v                  Effective Porosity
  SW          v/v                  Water Saturation
  SHC         v/v                  Hydrocarbon Saturation
  BVH         v/v                  Bulk Volume Hydrocarbon
  BVW         v/v                  Bulk Volume Water
  BVWI        v/v                  Bulk Volume Water Irreducible
  BVWF        v/v                  Bulk Volume Water Free
  BVOM        v/v                  Bulk Volume Fraction Organic Matter
  BVCLAY      v/v                  Bulk Volume Fraction Clay
  BVPYR       v/v                  Bulk Volume Fraction Pyrite
  VOM         v/v                  Matrix Volume Fraction Organic Matter
  VCLAY       v/v                  Matrix Volume Fraction Clay
  VPYR        v/v                  Matrix Volume Fraction Pyrite
  RHOM        g/cc                 Matrix Density
  TOC         wt/wt                Matrix Weight Fraction Organic Matter
  WTCLAY      wt/wt                Matrix Weight Fraction Clay
  WTPYR       wt/wt                Matrix Weight Fraction Pyrite
  BVQTZ       v/v                  Bulk Volume Fraction Quartz
  VQTZ        v/v                  Matrix Volume Fraction Quartz
  WTQTZ       wt/wt                Matrix Weight Fraction Quartz
  BVCLC       v/v                  Bulk Volume Fraction Calcite
  VCLC        v/v                  Matrix Volume Fraction Calcite
  WTCLC       wt/wt                Matrix Weight Fraction Calcite
  BVDOL       v/v                  Bulk Volume Fraction Dolomite
  VDOL        v/v                  Matrix Volume Fraction Dolomite
  WTDOL       wt/wt                Matrix Weight Fraction Dolomite
  OIP         wt/wt                Matrix Weight Fraction Dolomite
