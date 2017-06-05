# -*- coding: utf-8 -*-
"""
Log contains parent classes to work with log data, the Parameter class
and the Log class.

The Log class is subclassed based on the source of information. This way,
regardless of the data source, petrophysical calculations are available on the
logging data. The log data itself resides in a PANDAS dataframe for manipulation.

"""

import os
import re
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import datetime as dt
from scipy.optimize import nnls


class Parameter(object):
    """
    A Parameter from the log header data.

    The parameter class stores header data from logs in both version 1 and 2
    formats. There are two main ways to initialize:

        1. providing a header line from an las file

            This reads the line and parses out correct values for las file of
            version 1 and version 2

        2. providing specific paramter values

            This specifies each parameter value and leaves line as an empty
            str.

    Attributes
    ----------
    line: str
        Unformated raw str of the header line from an las file
    name: str
        name of the well parameter
    unit: str
        unit of the well parameter
    value: str
        value of the well parameter
    right_value: str
        right adjusted value for las version 1
    des: str
        parameter description
    """

    def __init__(self, line = '', name = '', unit = '', value = '', right_value = '', des = ''):

        """
        initialize Parameter object

        Parameters
        ----------
        line: str, optional
            str of the header line from an las file
        name: str, optional
            name of the well parameter
        unit: str, optional
            unit of the well parameter
        value: str, optional
            value of the well parameter
        right_value: str, optional
            right adjusted value for las version 1
        des: str, optional
            parameter description

        Raises
        ------
        ValueError
            If name or line is not specified, initialization raises a value
            error. The name or line variable is required to correctly parse
            into an object.

        Example
        -------
        >>> kb = Parameter(name = 'EKB', unit = 'FT', value = '1000.00', des = 'Elevation of Kelly Bushing')

        """

        if len(line) > 0:
            line = line.strip()
            self.name = line.split('.')[0].strip()
            start = len(line.split('.')[0]) + 1

            if len(line.split(':')) == 2:
                ### only one : in line ###
                self.des = line.split(':')[-1].strip()
            else:
                ### check for time stamp ###
                matches = re.compile(r'\d{2}:\d{2}').findall(line)
                if len(matches) > 0:
                    ### time stamp in value field ###
                    first_index = line.index(':')
                    second_index = line.index(':', first_index + 1)
                    self.des = line[second_index + 1:].strip()
                else:
                    ### : found in description ###
                    index = line.index(':') + 1
                    self.des = line[index:].strip()

            if len(self.des) > 3:
                end = len(line) - (len(self.des) + 3)
            else:
                end = len(line) - (len(self.des) + 1)

            line = line[start: end].split(' ')
            i = -1
            self.right_value = ''
            if line[-1] != '':
                right_values = []
                data = line[i]
                while data != '':
                    right_values.append(data)
                    i -= 1
                    data = line[i]
                for value in reversed(right_values):
                    self.right_value += value + ' '
                self.right_value = self.right_value[:-1]
            self.unit = line[0]
            line = list(filter(lambda x: x != '', line[1:i]))
            self.value = ''
            for l in line:
                self.value += l + ' '
            self.value = self.value[:-1].strip()
        elif len(name) > 0:
            self.name = name
            self.unit = unit
            self.value = value
            self.right_value = right_value
            self.des = des
        else:
            raise ValueError('A name or line must be specified to create a Parameter object.')

    def to_string(self):
        """
        Format parameter string for las file export.

        Returns
        ------
        output
            A string of the Parameter for use in las file format version 2

        Example
        -------
        >>> kb = Parameter(name = 'EKB', unit = 'FT', value = '1000.00', des = 'Elevation of Kelly Bushing')
        >>> print(kb.to_string())
        EKB         .F           1000.00                           :Elevation of Kelly Bushing

        """
        output = self.name.ljust(12, ' ') + '.'
        output += self.unit.ljust(12, ' ')
        output += self.value.ljust(36, ' ')
        right_adjustment = max(0 - (len(self.right_value) + 2), -44)
        if right_adjustment != -2:
            output = output[:right_adjustment] + self.right_value + ': '
        else:
            output = output[:right_adjustment] + self.right_value + ':'
        output += self.des
        return output


class Log(object):
    """
    Well Log

    Parent class for all petrophysical calculations with subclasses for reading
    and writing to specific databases or data types. Parameter names
    are stored in a list to maintain order as order is important when reading
    and writing curve data.

    Attributes
    ----------
    uwi : str
        Unique Well Identifier
    null : str
        null value in log
    version_parameters : list
        A list of the names of parameters assocaited with the version of the las
        file.
    version_values : dict
        A dictionary of Parameter objects associated with the las version
        section of the las file. The names of the objects are keys and found in
        the version_parameters list.
    well_parameters : list
        A list of the names of parameters assocaited with the well in the las
        header.
    well_values : dict
        A dictionary of Parameter objects associated with the well parameter
        section of the las file. The names of the objects are keys and found in
        the well_parameters list.
    parameter_parameters : list
        A list of the names of parameters assocaited with the parameter section
        in the las header.
    parameter_values : dict
        A dictionary of Parameter objects associated with the parameters
        section of the las file. The names of the objects are keys and found in
        the parameter_parameters list.
    curve_parameters : list
        A list of the names of curves in the las file.
    curve_values : dict
        A dictionary of Parameter objects associated with the curve
        section of the las file. The names of the objects are keys and found in
        the curve_parameters list.
    curve_df : DataFrame
        A PANDAS DataFrame containing the log responses
    fluid_properties_parameters : dict
        A dictionary of parameter dictionaries for fluid properties calculations
    multimineral_parameters : dict
        A dictionary of parameter dictionaries for multimineral calculations
    tops : dict
        A dictionary of tops : depth


    See also
    --------
    Las
        The subclass of Log for working with raw or standalone las files.

    """


    def __init__(self):

        self.uwi = ''
        self.null = ''

        self.version_parameters = []
        self.version_values = {}

        self.well_parameters = []
        self.well_values = {}

        self.parameter_parameters = []
        self.parameter_values = {}

        self.curve_parameters = []
        self.curve_values = {}

        self.curve_df = pd.DataFrame()

        self.fluid_properties_parameters = {}
        self.multimineral_parameters = {}

        self.tops = {}


    def precondition(self, drho_matrix = 2.71):
        """
        Preconditions log curve by aliasing names.

        Precondition is used after initializing data by subclasses and
        standardizes names for future calculations.

        Notes
        -----
            1. Curve Alias is provided by the curve_alias.xml file


        See Also
        --------
        Las : :__init__: ~`las.Las.__init__`
            Uses precondition by default when reading las file

        Parameters
        ----------
        drho_matrix : float, optional
            drho_matrix is for converting density porosity to bulk densty, and
            is only used when bulk density is missing. Default value for
            limestone matrix. If log was run on sandstone matrix, use 2.65. If
            log was run on dolomite matrix, use 2.85.

        """

        file_dir = os.path.dirname(__file__)
        ALIAS_XML_PATH = os.path.join(file_dir, '..', 'data', 'sample', 'curve_alias.xml')

        if not hasattr(self, 'curve_df'):
            raise ValueError('curve_df not found for log')

        if not hasattr(self, 'curve_parameters'):
            raise ValueError('Curve Names not found for log')

        if not os.path.isfile(ALIAS_XML_PATH):
            raise ValueError('Could not find alias xml at: %s' % ALIAS_XML_PATH)

        with open(ALIAS_XML_PATH, 'r') as f:
            root = ET.fromstring(f.read())

        for alias in root:
            for curve in alias:
                if curve.tag in self.curve_parameters:
                    if alias.tag not in self.curve_parameters:
                        self.curve_df[alias.tag] = self.curve_df[curve.tag]
                        line = self.curve_values[curve.tag]
                        name = alias.tag
                        unit = line.unit
                        des = line.des
                        value = line.value
                        right_value = line.right_value
                        self.curve_parameters.append(alias.tag)
                        self.curve_values[name] = Parameter(name = name,
                                                            unit = unit,
                                                            value = value,
                                                            right_value = right_value,
                                                            des = des)
                    break

        if 'RHOB_N' not in self.curve_parameters and 'DPHI_N' in self.curve_parameters:
            self.curve_parameters.append('RHOB_N')
            self.curve_values['RHOB_N'] = Parameter(name = 'RHOB_N',
                                                    unit = 'g/cc',
                                                    des = 'Bulk Density from DPHI')
            self.curve_df['RHOB_N'] = drho_matrix - (drho_matrix - 1) * self.curve_df.DHI_N

        ### rename depth track to DEPTH for consistency ###
        depth_name = self.curve_parameters[0]
        self.curve_values[depth_name].name = 'DEPTH'
        self.curve_df.columns.values[0] = 'DEPTH'

        self.curve_df.replace(self.null, np.nan, inplace = True)


    def tops_from_csv(self, csv_path = None):
        """
        Reads tops from a csv file and saves as dictionary.

        Parameters
        ----------
        csv_path : str (default None)
            Path to csv file to read. Must contain header row the the following
            properties:
                uwi : str
                    Unique Well Identifier
                form : str
                    Name of formation top
                depth : float
                    depth of corresponding formation top

        Notes
        -----

        Format for csv:

            uwi,form,depth
            11111111111,WFMPA,7000.50
            11111111111,WFMPB,7250.50
            11111111111,WFMPC,7500.00
            11111111111,WFMPD,7700.25
            11111111111,DEAN,8000.00

        Examples
        --------
        >>> import petropy as ptr
        >>> log = ptr.logdataset('WFMP') # loads Wolfcamp Log
        >>> log.tops_from_csv() # loads default tops for included datasets

        >>> import petropy as ptr
        >>> log = ptr.Las('path/to/well.las') # loads specified las file
        >>> log.tops_from_csv('path/to/tops.csv') # loads specified tops csv

        """

        if csv_path is None:
            local_path = os.path.dirname(__file__)
            csv_path = os.path.join(local_path, '..', 'data', 'sample', 'tops.csv')

        top_df = pd.read_csv(csv_path, dtype = {'uwi': str, 'form': str, 'depth': float})

        for r, row in top_df[top_df.uwi == self.uwi].iterrows():
            self.tops[row.form] = row.depth


    def next_formation_depth(self, formation):
        """
        Return top of formation below specified formation.

        Function required since python 2 does not use ordered dictionaries.

        Parameter
        --------
        formation : str
            name of formation which should be found in preloaded tops

        Returns
        -------
        bottom : float
            top of formation below input formation, acting as bottom of input
            formation

        Example
        -------
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> log.tops_from_csv() # loads default tops for included datasets
        >>> wfmpa_top = log.tops['WFMPA']
        >>> print(wfmpa_top)
        6993.5
        >>> wfmpa_bottom = log.next_formation_depth('WFMPB')
        >>> print(wfmpa_bottom)
        7294.0
        >>> wfmpb_top = log.tops['WFMPB']
        >>> print(wfmpb_top)
        7294.0

        """

        top = self.tops[formation]

        bottom = self.curve_df.DEPTH.max()
        closest_formation = bottom - top
        for form in self.tops:
            form_depth = self.tops[form]
            if form_depth > top and form_depth - top < closest_formation:
                bottom = form_depth
                closest_formation = bottom - top

        return bottom


    def fluid_properties_parameters_from_csv(self, csv_path = None):
        """
        Reads parameters from a csv for input into fluid properties.

        This method reads the file located at the csv_path and turns the values
        into dictionaries to be used as inputs into the fluid_properties method.

        Parameters
        ----------
        csv_path : str (default None)
            Path to csv file to read. Must contain header row the the following
            properties:
                mast : float (default 67)
                    The mean annual surface temperature at the location of the
                    well in degrees Fahrenheit.
                temp_grad : float (default 0.015)
                    The temperature gradient of the reservoir in °F / ft.
                press_grad : float (default 0.5)
                    The pressure gradient of the reservoir in psi / ft.
                rws : float (default 0.1)
                    The resistivity of water at surface conditions in ohm.m.
                rwt : float (default 70)
                    The temperature of the rws measurement in °F.
                rmfs : float (default 0.4)
                    The resistivity of mud fultrate at surface conditions in
                    ohm.m.
                rmft : float (default 100)
                    The temperature of the rmfs measurement in °F.
                gas_grav : float (default 0.67)
                    The specific gravity of the separator gas. Air = 1,
                    CH4 = 0.577
                oil_api : float (default 38)
                    The api gravity of oil after the separator. If fluid system
                    is dry gas, use oil_api = 0.
                p_sep : float (default 100)
                    The pressure of the separator, assuming a 2 stage system.
                    Only used when oil_api is > 0 (not dry gas).
                t_sep : float
                    The temperature of the separator , assuming a 2 stage
                    system. Only used with oil_api > 0.
                yn2 : float (default 0)
                    Molar fraction of nitrogren in gas.
                yco2 : float (default 0)
                    Molar fration of carbon dioxide in gas.
                yh2s : float (default 0)
                    Molar fraction of hydrogen sulfide in gas.
                yh20: float (default 0)
                    Molar fraction of water in gas.
                rs : float (default 0)
                    Solution gas oil ratio at reservoir conditions. If unknwon,
                    use 0 and correlation will be calculated.
                lith_grad: float (default 1.03)
                    Lithostatic gradient in psi / ft.
                biot : float (default 0.8)
                    Biot constant.
                pr : float (default 0.25)
                    Poissons ratio

        Examples
        --------
        >>> import petropy as ptr
        >>> from petropy import datasets
        >>> log = datasets('WFMP') # loads Wolfcamp Log
        >>> log.fluid_properties_parameters_from_csv() # loads base parameters

        >>> import petropy as ptr
        >>> from petropy import datasets
        >>> log = datasets('WFMP') # loads Wolfcamp Log
        >>> my_csv_paramters = 'path/to/csv/file.csv' # specified csv
        >>> log.fluid_properties_parameters_from_csv(my_csv_paramters) # loads specified parameters

        See Also
        --------
        fluid_properties
            Calculates fluid properties using input settings loaded through this
            method

        """

        if csv_path is None:
            local_path = os.path.dirname(__file__)
            csv_path = os.path.join(local_path, '..', 'data', 'sample', 'fluid_properties_parameters.csv')

        param_df = pd.read_csv(csv_path)
        param_df = param_df.set_index('name')

        self.fluid_properties_parameters = param_df.to_dict(orient = 'index')


    def fluid_properties(self, top = None, bottom = None, mast = 67, temp_grad = 0.015, press_grad = 0.5, rws = 0.1, rwt = 70, rmfs = 0.4, rmft = 100, gas_grav = 0.67, oil_api = 38, p_sep = 100, t_sep = 100, yn2 = 0, yco2 = 0, yh2s = 0, yh20 = 0, rs = 0, lith_grad = 1.03, biot = 0.8, pr = 0.25):
        """
        Calculates fluid properties along wellbore.

        Current single phase fluid properties assumes either:

            1. Dry Gas at Reservoir Conditions
            Methane as hydrocarbon type with options to include N2, CO2, H2S, or
            H2O. To assume dry_gas, set oil_api = 0.

            2. Oil at Reservoir Conditions
            Assumes reservoir fluids are either a black oil or volatile oil.
            Separator conditions of gas are used to calculate bubble point and
            the reservoir fluid properties of the reconstituted oil.

        The output updates the curve_df attribute of the Log with the following
        calculated curves at each depth:

        PORE_PRESS : (psi)
            Reservoir pore pressure
        RES_TEMP : (°F)
            Reservoir temperature
        NES : (psi)
            Reservoir net effective stress
        RW : (ohm.m)
            Resistivity of water
        RMF : (ohm.m)
            Resistivity of mud filtrate
        RHO_HC : (g / cc)
            Density of hydrocarbon
        RHO_W : (g / cc)
            Density of formation water
        RHO_MF : (g / cc)
            Density of mud filtrate
        NPHI_HC
            Neutron log response of hydrocarbon
        NPHI_W
            Neutron log response of water
        NPHI_MF
            Neutron log response of mud filtrate
        MU_HC : (cP)
            Viscosity of hydrocarbon
        Z
            Compressiblity factor for non-ideal gas. Only output if oil_api = 0
        CG : (1 / psi)
            Gas Compressiblity. Only output if oil_api = 0
        BG
            Gas formation volume factor. Only output if oil_api = 0
        BP : (psi)
            Bubble point. Only output if oil_api > 0
        BO
            Oil formation volume factor. Only output if oil_api > 0

        Parameters
        ----------
        top : float, optional
            The top depth to begin fluid properties calculation. If value is not
            specified, the calculations will start at the top of the log.
        bottom : float, optional
            The bottom depth to end fluid properties, inclusive. If the value is
            not specified, the calcuations will go to the end of the log.
        mast : float (default 67)
            The mean annual surface temperature at the location of the well in
            degrees Fahrenheit.
        temp_grad : float (default 0.015)
            The temperature gradient of the reservoir in °F / ft.
        press_grad : float (default 0.5)
            The pressure gradient of the reservoir in psi / ft.
        rws : float (default 0.1)
            The resistivity of water at surface conditions in ohm.m.
        rwt : float (default 70)
            The temperature of the rws measurement in °F.
        rmfs : float (default 0.4)
            The resistivity of mud fultrate at surface conditions in ohm.m.
        rmft : float (default 100)
            The temperature of the rmfs measurement in °F.
        gas_grav : float (default 0.67)
            The specific gravity of the separator gas. Air = 1, CH4 = 0.577
        oil_api : float (default 38)
            The api gravity of oil after the separator. If fluid system is dry
            gas, use oil_api = 0.
        p_sep : float (default 100)
            The pressure of the separator, assuming a 2 stage system. Only used
            when oil_api is > 0 (not dry gas).
        t_sep : float
            The temperature of the separator , assuming a 2 stage system. Only
            used with oil_api > 0.
        yn2 : float (default 0)
            Molar fraction of nitrogren in gas.
        yco2 : float (default 0)
            Molar fration of carbon dioxide in gas.
        yh2s : float (default 0)
            Molar fraction of hydrogen sulfide in gas.
        yh20: float (default 0)
            Molar fraction of water in gas.
        rs : float (default 0)
            Solution gas oil ratio at reservoir conditions. If unknwon, use 0
            and correlation will be used.
        lith_grad: float (default 1.03)
            Lithostatic overburden pressure gradient in psi / ft.
        biot : float (default 0.8)
            Biot constant.
        pr : float (default 0.25)
            Poissons ratio


        References
        ----------
        Ahmed, Tarek H. Reservoir Engineering Handbook. Oxford: Gulf
            Professional, 2006.

        Lee, John, and Robert A. Wattenbarger. Gas Reservoir Engineering.
            Richardson, TX: Henry L. Doherty Memorial Fund of AIME, Society of
            Petroleum Engineers, 2008.


        Examples
        --------
        >>> import petropy as ptr
        >>> from petropy import datasets
        >>> log = datasets('WFMP') # loads Wolfcamp Log
        >>> log.fluid_properties() # calculates fluid properties with default settings


        See Also
        --------
        fluid_properties_parameters_from_csv
            loads properties from preconfigured csv file
        multimineral_model
            builds on fluid properties to calculate full petrophysical model

        """

        if 'DEPTH' not in list(self.curve_df.columns.values):
            raise ValueError('No DEPTH curve found in log.')

        if top is not None:
            df = self.curve_df[self.curve_df.DEPTH >= top].copy()
        else:
            df = self.curve_df.copy()
            top = self.curve_df.DEPTH.min()

        if bottom is not None:
            df = df[df.DEPTH <= bottom].copy()
        else:
            bottom = self.curve_df.DEPTH.max()

        ### fluid property calculations ###

        form_temp = mast + temp_grad * df.DEPTH
        pore_press = press_grad * df.DEPTH

        ### water properties ###
        rw = (rwt + 6.77) / (form_temp + 6.77) * rws
        rmf = (rmft + 6.77) / (form_temp + 6.77) * rmfs

        rw68 = (rwt + 6.77) / (68 + 6.77) * rws
        rmf68 = (rmft + 6.77) / (68 + 6.77) * rws

        ### weight percent total disolved solids ###
        xsaltw = 10 ** (-0.5268 * (np.log10(rw68) ) ** 3 - 1.0199 * (np.log10(rw68)) ** 2 - 1.6693 * (np.log10(rw68)) - 0.3087)
        xsaltmf = 10 ** (-0.5268 * (np.log10(rmf68) ) ** 3 - 1.0199 * (np.log10(rmf68)) ** 2 - 1.6693 * (np.log10(rmf68)) - 0.3087)

        ### bw for reservoir water. Eq 1.83 - 1.85 Gas Reservoir Engineering ###
        dvwt = -1.0001 * 10 ** -2 + 1.33391 * 10 ** -4 * form_temp + 5.50654 * 10 ** -7 * form_temp ** 2
        dvwp = -1.95301 * 10 ** -9 * pore_press * form_temp - 1.72834 * 10 ** -13 * pore_press ** 2 * form_temp - 3.58922 * 10 ** -7 * pore_press - 2.25341 * 10 ** -10 * pore_press ** 2
        bw = (1 + dvwt) * (1 + dvwp)

        ### calculate solution gas in water ratio, Eq. 1.86 - 1.91 Gas Reservoir Engineering ###
        rsa = 8.15839 - 6.12265 * 10 ** -2 * form_temp + 1.91663 * 10 **-4 * form_temp ** 2 - 2.1654 * 10 ** -7 * form_temp ** 3
        rsb = 1.01021 * 10 ** -2 - 7.44241 * 10 ** -5 * form_temp + 3.05553 * 10 ** -7 * form_temp ** 2 - 2.94883 * 10 ** -10 * form_temp ** 3
        rsc = -1.0 * 10 ** -7 * (9.02505 - 0.130237 * form_temp + 8.53425 * 10 ** -4 * form_temp ** 2 - 2.34122 * 10 ** -6 * form_temp ** 3 + 2.37049 * 10 ** -9 * form_temp ** 4)
        rswp = rsa + rsb * pore_press + rsc * pore_press ** 2
        rsw = rswp * 10 ** (-0.0840655 * xsaltw * form_temp ** -0.285584)

        ### log responses ###
        rho_w = (2.7512 * 10 ** -5 * xsaltw + 6.9159 * 10 ** -3 * xsaltw + 1.0005) * bw
        rho_mf = (2.7512 * 10 ** -5 * xsaltmf + 6.9159 * 10 ** -3 * xsaltmf + 1.0005) * bw
        nphi_w = 1 + 0.4 * (xsaltw / 100)
        nphi_mf = 1 + 0.4 * (xsaltmf / 100)

        ### net efective stress ###
        nes = (((lith_grad * df.DEPTH) - (biot * press_grad * df.DEPTH) + 2 * (pr / (1 - pr)) * (lith_grad * df.DEPTH) - (biot * press_grad * df.DEPTH))) / 3

        ### gas reservoir ###
        if oil_api == 0:
            # hydrocarbon garvity only
            hc_grav = (gas_grav - 1.1767 * yh2s - 1.5196 * yco2 - 0.9672 * yn2 - 0.622 * yh20) / (1.0 - yn2 - yco2 - yh20 - yh2s)
            # pseudocritical properties of hydrocarbon
            ppc_h = 756.8 - 131.0 * hc_grav - 3.6 * (hc_grav ** 2)
            tpc_h = 169.2 + 349.5 * hc_grav - 74.0 * (hc_grav ** 2)
            # pseudocritical properties of mixture
            ppc = (1.0 - yh2s - yco2 - yn2 - yh20) * ppc_h + 1306.0 * yh2s + 1071.0 * yco2 + 493.1 * yn2 + 3200.1 * yh20
            tpc = (1.0 - yh2s - yco2 - yn2 - yh20) * tpc_h + 672.35 * yh2s + 547.58 * yco2 + 227.16 * yn2 + 1164.9 * yh20

            # Wichert-Aziz correction for hydrogen sulfide and carbon dioxide
            if yco2 > 0 or yh2s > 0:
                epsilon = 120 * ((yco2 + yh2s) ** 0.9 - (yco2 + yh2s) ** 1.6) + 15 * (yh2s ** 0.5 - yh2s ** 4)
                tpc_temp = tpc - epsilon
                ppc = (ppc_a * tpc_temp) / (tpc + (yh2s * (1.0 - yh2s) * epsilon))
                tpc = tpc_temp
            # Casey correction for nitrogen and water vapor
            if yn2 > 0 or yh20 > 0:
                tpc_cor = -246.1 * yn2 + 400 * yh20
                ppc_cor = -162.0 * yn2 + 1270.0 * yh20
                tpc = (tpc - 227.2 * yn2 - 1165.0 * yh20) / (1.0 - yn2 - yh20) + tpc_cor
                ppc = (ppc - 493.1 * yn2 - 3200.0 * yH20) / (1.0 - yn2 - yh20) + ppc_cor
            # Reduced pseudocritical properties
            tpr = (form_temp + 459.67) / tpc
            ppr = pore_press / ppc

            ### z factor from Dranchuk and Abou-Kassem fit of Standing and Katz chart ###
            a = [0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.721]

            t2 = a[0] * tpr + a[1] + a[2] / (tpr ** 2) + a[3] / (tpr ** 3) + a[4] / (tpr ** 4)
            t3 = a[5] * tpr + a[6] + a[7] / tpr
            t4 = -a[8] * (a[6] + a[7] / tpr)
            t5 = a[9] / (tpr ** 2)
            r = 0.27 * ppr / tpr
            counter = 0
            diff = 1
            z = 0.27 * ppr / tpr / r
            while counter <= 10 and diff > 10 ** -5:
                counter += 1
                f = r * (tpr + t2 * r + t3 * r ** 2 + t4 * r ** 5 + t5 * r ** 2 * (1 + a[10] * r**2) * np.exp(-a[10] * r **2)) - 0.27 * ppr
                fp = tpr + 2 * t2 * r + 3 * t3 * r ** 2 + 6 * t4 * r ** 5 + t5 * r ** 2 * np.exp(-a[10] * r ** 2) * (3 + a[10] * r ** 2 * (3 - 2 * a[10] * r ** 2))
                r = r - f/fp
                diff = np.abs(z - (0.27 * ppr / tpr / r)).max()
                z = 0.27 * ppr / tpr / r

            ### gas compressiblity from Dranchuk and Abau-Kassem ###
            cpr = tpr * z / ppr / fp
            cg = cpr / ppc

            ### gas expansion factor ###
            bg = (0.0282793 * z * (form_temp + 459.67)) / pore_press

            ### gas density Eq 1.64 GRE ###
            rho_hc = 1.495 * 10 ** -3 * (pore_press * (gas_grav)) / (z * (form_temp + 459.67))
            nphi_hc = 2.17 * rho_hc

            ### gas viscosity Lee Gonzalez Eakin method Eqs. 1.63-1.67 GRE ###
            k = ((9.379 + 0.01607 * (28.9625 * gas_grav)) * (form_temp + 459.67) ** 1.5) / (209.2 + 19.26 * (28.9625 * gas_grav) + (form_temp + 459.67))
            x = 3.448 + 986.4 / (form_temp + 459.67) + 0.01009 * (28.9625 * gas_grav)
            y = 2.447 - 0.2224 * x
            mu_hc = 10 **-4 * k * np.exp(x * rho_hc ** y)


        ### oil reservoir ###
        else:

            # Normalize gas gravity to separator pressure of 100 psi
            ygs100 = gas_grav * (1 + 5.912 * 0.00001 * oil_api * (t_sep - 459.67) * np.log10(p_sep / 114.7))

            if oil_api < 30:
                if rs == 0 or rs is None:
                    rs = 0.0362 * ygs100 * pore_press ** 1.0937 * np.exp((25.724 * oil_api) / (form_temp + 459.67))
                bp = ((56.18 * rs / ygs100) * 10 ** (-10.393 * oil_api / (form_temp + 459.67))) ** 0.84246
                ### gas saturated bubble-point ###
                bo = 1 + 4.677 * 10 ** -4 * rs + 1.751 * 10 ** -5 * (form_temp - 60) * (oil_api / ygs100) - 1.811 * 10 ** -8 * rs * (form_temp - 60) * (oil_api / ygs100)
            else:
                if rs == 0 or rs is None:
                    rs = 0.0178 * ygs100 * pore_press ** 1.187 * np.exp((23.931 * oil_api) / (form_temp + 459.67))
                bp = ((56.18 * rs / ygs100) * 10 ** (-10.393 * oil_api / (form_temp + 459.67))) ** 0.84246
                ### gas saturated bubble-point ###
                bo = 1 + 4.670 * 10 ** -4 * rs + 1.1 * 10 ** -5 * (form_temp - 60) * (oil_api / ygs100) + 1.337 * 10 ** -9 * rs * (form_temp - 60) * (oil_api / ygs100)

            ### calculate bo for undersaturated oil ###
            pp_bp = pore_press > bp + 100
            bo[pp_bp] = bo[pp_bp] * np.exp(-(0.00001 * (-1433 + 5 * rs + 17.2 * form_temp[pp_bp] - 1180 * ygs100 + 12.61 * oil_api)) * np.log(pore_press[pp_bp] / bp[pp_bp]))

            ### oil properties ###
            rho_hc = (((141.5 / (oil_api + 131.5) * 62.428) + 0.0136 * rs *ygs100) / bo) / 62.428
            nphi_hc = 1.003 * rho_hc

            ### oil viscosity from Beggs-Robinson, RE Handbook Eqs. 2.121
            muod = 10 ** (np.exp(6.9824 - 0.04658 * oil_api) * form_temp ** -1.163) - 1
            mu_hc = (10.715 * (rs + 100) ** -0.515) * muod ** (5.44 * (rs + 150) ** -0.338)

            # undersaturated oil viscosity, Vasquez and Beggs Eqs. 2.123
            mu_hc[pp_bp] = mu_hc[pp_bp] * (pore_press[pp_bp] / bp[pp_bp]) ** (2.6 * pore_press[pp_bp] ** 1.187 * 10 ** (-0.000039 * pore_press[pp_bp] - 5))


        ### format outputs into curve_df ###

        output_curves = {
            'PORE_PRESS': pore_press,
            'RES_TEMP': form_temp,
            'NES': nes,
            'RW': rw,
            'RMF': rmf,
            'RHO_HC': rho_hc,
            'RHO_W': rho_w,
            'RHO_MF': rho_mf,
            'NPHI_HC': nphi_hc,
            'NPHI_W': nphi_w,
            'NPHI_MF': nphi_mf,
            'MU_HC': mu_hc,
        }

        output_curve_parameters = {
            'PORE_PRESS': Parameter(name = 'PORE_PRESS', unit = 'psi', des = 'Calculated Pore Pressure'),
            'RES_TEMP': Parameter(name = 'RES_TEMP', unit = 'F', des = 'Calculated Reservoir Temperature'),
            'NES': Parameter(name = 'NES', unit = 'psi', des = 'Calculated Net Effective Stress'),
            'RW': Parameter(name = 'RW', unit = 'ohmm', des = 'Calculated Resistivity Water'),
            'RMF': Parameter(name = 'RMF', unit = 'ohmm', des = 'Calculated Resistivity Mud Filtrate'),
            'RHO_HC': Parameter(name = 'RHO_HC', unit = 'g/cc', des = 'Calculated Density of Hydrocarbon'),
            'RHO_W': Parameter(name = 'RHO_W', unit = 'g/cc', des = 'Calculated Density of Water'),
            'RHO_MF': Parameter(name = 'RHO_MF', unit = 'g/cc', des = 'Calculated Density of Mud Filtrate'),
            'NPHI_HC': Parameter(name = 'NPHI_HC', unit = 'v/v', des = 'Calculated Neutron Log Response of Hydrocarbon'),
            'NPHI_W': Parameter(name = 'NPHI_W', unit = 'v/v', des = 'Calculated Neutron Log Response of Water'),
            'NPHI_MF': Parameter(name = 'NPHI_MF', unit = 'v/v', des = 'Calculated Neutron Log Response of Mud Filtrate'),
            'MU_HC': Parameter(name = 'MU_HC', unit = 'cP', des = 'Calculated Viscosity of Hydrocarbon')
        }

        for curve in output_curves:
            self.curve_parameters.append(curve)
            self.curve_values[curve] = output_curve_parameters[curve]
            self.curve_df.loc[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom), curve] = output_curves[curve]

        ### gas curves ###
        if oil_api == 0:
            gas_curves = {
                'Z': z,
                'CG': cg,
                'BG': bg
            }

            gas_curve_parameters = {
                'Z': Parameter(name = 'Z', des = 'Calculated Real Gas Z Factor'),
                'CG': Parameter(name = 'CG', unit = '1 / psi', des = 'Calculated Gas compressibility'),
                'BG': Parameter(name = 'BG', des = 'Calculated Gas Formation Volume Factor')
            }

            for curve in gas_curves:
                self.curve_parameters.append(curve)
                self.curve_values[curve] = gas_curve_parameters[curve]
                self.curve_df.loc[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom), curve] = gas_curves[curve]

        ### oil curves ###
        else:
            oil_curves = {
                'BO': bo,
                'BP': bp
            }

            oil_curve_parameters = {
                'BO': Parameter(name = 'BO', des = 'Calculated Oil Formation Volume Factor'),
                'BP': Parameter(name = 'BP', unit = 'psi', des = 'Calculated Bubble Point')
            }

            for curve in oil_curves:
                self.curve_parameters.append(curve)
                self.curve_values[curve] = oil_curve_parameters[curve]
                self.curve_df.loc[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom), curve] = oil_curves[curve]


    def formation_fluid_properties(self, formations = [], parameter = 'default'):
        """
        Calculate fluid properties over formations with loaded paramters

        Parameters
        ----------
        formations : list (default [])
            list of formations to calculate fluid properties over
        parameter : str (default 'default')
            name of parameter to use for fluid properties parameter settings
            loaded in method fluid_properties_parameters_from_csv

        See Also
        --------
        fluid_properties
            calculates fluid properties of log
        fluid_properties_parameters_from_csv
            loads fluid properties parameters
        tops_from_csv
            loads tops of log from csv

        """

        for form in formations:
            top = self.tops[form]
            bottom = self.next_formation_depth(form)

            params = self.fluid_properties_parameters[parameter]

            self.fluid_properties(top = top, bottom = bottom, **params)


    def multimineral_parameters_from_csv(self, csv_path = None):
        """
        Reads parameters from a csv for input into the multimineral model.

        This method reads the file located at the csv_path and turns the values
        into dictionaries to be used as inputs into the multimineral method.

        Parameters
        ----------
        csv_path : str (default None)
            Path to csv file to read. Must contain header row with the following
            properties:

            gr_matrix : float (default 10)
                Gamma Ray response of clean (non-clay) matrix
            nphi_matrix : float (default 0)
                Neutron response of clean (non-clay) matrix
            gr_clay : float (default 450)
                Gamma Ray response of pure clay matrix
            rho_clay : float (default 2.64)
                Density of pure clay matrix
            nphi_clay : float (default 0.65)
                Neutron reponse of pure clay matrix
            pe_clay : float (default 4)
                Photoelectric response of pure clay matrix
            rma : float (default 180)
                Resistivity of clean tight matrix
            rt_clay : float (default 80)
                Resistivity for inorganic shale matrix
            vclay_linear_weight : float (default 1)
                Weight of liner clay volume
            vclay_clavier_weight : float (defaul 0.5)
                Weight of Clavier clay volume
            vclay_larionov_weight : float (defaul 0.5)
                Weight of Larionov clay volume
            vclay_nphi_weight : float (default 1)
                Weight of Neutron clay volume
            vclay_nphi_rhob_weight : float (default 1)
                Weight of Neutron Density clay volume
            vclay_cutoff : float (default 0.05)
                Cutoff for oranics calculation. If vclay < vclay_cutoff then toc = 0
            rho_om : float (default 1.15)
                Density of organic matter
            nphi_om : float (default 0.6)
                Neutron response of pure organic matter
            pe_om : float (default 0.2)
                Photoelectric response of pure organic matter
            ro : float (default 1.6)
                Vitronite reflectance of organic matter
            lang_press : float (default 670)
                Langmiur pressure for gas adsorption on organic matter in psi
            passey_nphi_weight : float (default 1)
                Weight for Passey nphi toc
            passey_rhob_weight : float (default 1)
                Weight for Passey rhob toc
            passey_lom : float (default 10)
                Passey level of organic maturity
            passey_baseline_res : float (default 40)
                Passey inorganic baseline resistivity
            passey_baseline_rhob : float (default 2.65)
                Passey inorganic baseline density
            passey_baseline_nphi : float (default 0)
                Passey inorganic baseline neutron
            schmoker_weight : float (default 1)
                Weight for Schmoker toc
            schmoker_slope : float (default 0.7257)
                Slope for schmoker density to toc correlation
            schmoker_baseline_rhob : float (default 2.6)
                Density cutoff for schmoker toc correlation
            rho_pyr : float (default 5)
                Density of pyrite
            nphi_pyr : float (default 0.13)
                Neutron response of pure pyrite
            pe_pyr : float (default 13)
                Photoelectric response of pure pyrite
            om_pyrite_slope : float (default 0.2)
                Slope correlating pyrite volume to organic matter
            include_qtz : str {'YES', 'NO'} (default 'YES')
                Toggle to include or exclude qtz. 'YES' to include. 'NO' to exclude.
            rho_qtz : float (default 2.65)
                Density of quartz
            nphi_qtz : float (default -0.04)
                Neutron response for pure quartz
            pe_qtz : float (default 1.81)
                Photoelectric response for pure quartz
            include_clc : str {'YES', 'NO'} (default 'YES')
                Toggle to include or exclude clc. 'YES' to include. 'NO' to exclude.
            rho_clc : float (default 2.71)
                Density of calcite
            nphi_clc : float (default 0)
                Neutron response for pure calcite
            pe_clc : float (default 5.08)
                Photoelectric response for pure calcite
            include_dol : str {'YES', 'NO'} (default 'YES')
                Toggle to include or exclude dol. 'YES' to include. 'NO' to exclude.
            rho_dol : float (default 2.85)
                Density of dolomite
            nphi_dol : float (default 0.04)
                Neutron response to dolomite
            pe_dol : float (default 3.14)
                Photoelectric response to dolomite
            include_x : str {'YES', 'NO'} (default 'NO')
                Toggle to include or exclude exotic mineral, x. 'YES' to include.
                'NO' to exclude.
            name_x : str (default 'Gypsum')
                Name of exotic mineral, x.
            name_log_x : str (default 'GYP')
                Log name of exotic mineral, x
            rho_x : float (default 2.35)
                Density of exotic mineral, x
            nphi_x : float (default 0.507)
                Neutron response of exotic mineral, x
            pe_x : float (default 4.04)
                Photoelectric respone of exotic mineral, x
            pe_fl : float (default 0)
                Photoelectric response of reservoir fluid.
            m : float (default 2)
                Cementation exponent
            n : float (default 2)
                Saturation exponent
            a : float (default 1)
                Cementation constant
            cec : float (default -1)
                Cation Exchange Capaticy for use in Waxman Smits Sw equation. If
                cec = -1, correlation equation is used to calculate cec.
            archie_weight : float (default 0)
                Weight for archie Sw
            indonesia_weight : float (default 1)
                Weight for Indonesia Sw
            simandoux_weight : float (default 0)
                Weight for Simandoux Sw
            modified_simandoux_weight : float (default 0)
                Weight for Modified Simandoux Sw
            waxman_smits_weight : float (default 0)
                Weight for Waxman Smits Sw
            buckles_parameter : float (default -1)
                Buckles parameter for calculating irreducible water saturation. If
                less than 0, it is calculated using a correlation.

        Examples
        --------
        >>> import petropy as ptr
        >>> from petropy import datasets
        >>> log = datasets('WFMP') # loads Wolfcamp Log
        >>> log.multimineral_parameters_from_csv() # loads base parameters

        >>> import petropy as ptr
        >>> from petropy import datasets
        >>> log = datasets('WFMP') # loads Wolfcamp Log
        >>> my_csv_paramters = 'path/to/csv/file.csv' # specified csv
        >>> log.multimineral_parameters_from_csv(my_csv_paramters) # loads specified parameters

        See Also
        --------
        fluid_properties
            Calculates fluid properties using input settings loaded through this
            method

        """

        if csv_path is None:
            local_path = os.path.dirname(__file__)
            csv_path = os.path.join(local_path, '..', 'data', 'sample', 'multimineral_parameters.csv')

        param_df = pd.read_csv(csv_path)
        param_df = param_df.set_index('name')

        self.multimineral_parameters = param_df.to_dict(orient = 'index')


    def multimineral_model(self, top = None, bottom = None, gr_matrix = 10, nphi_matrix = 0, gr_clay = 350, rho_clay = 2.64, nphi_clay = 0.65, pe_clay = 4, rma = 180, rt_clay = 80, vclay_linear_weight = 1, vclay_clavier_weight = 0.5, vclay_larionov_weight = 0.5, vclay_nphi_weight = 1, vclay_nphi_rhob_weight = 1, vclay_cutoff = 0.1, rho_om = 1.15, nphi_om = 0.6, pe_om = 0.2, ro = 1.6, lang_press = 670, passey_nphi_weight = 1, passey_rhob_weight = 1, passey_lom = 10, passey_baseline_res = 40, passey_baseline_rhob = 2.65, passey_baseline_nphi = 0, schmoker_weight = 1, schmoker_slope =  0.7257, schmoker_baseline_rhob = 2.6, rho_pyr = 5, nphi_pyr = 0.13, pe_pyr = 13, om_pyrite_slope = 0.2, include_qtz = 'YES', rho_qtz = 2.65, nphi_qtz = -0.04, pe_qtz = 1.81, include_clc = 'YES', rho_clc = 1.71, nphi_clc = 0, pe_clc = 5.08, include_dol = 'YES', rho_dol = 2.85, nphi_dol = 0.04, pe_dol = 3.14, include_x = 'NO', name_x = 'Gypsum', name_log_x = 'GYP', rho_x = 2.35, nphi_x = 0.507, pe_x = 4.04, pe_fl = 0, m = 2, n = 2, a = 1, archie_weight = 0, indonesia_weight = 1, simandoux_weight = 0, modified_simandoux_weight = 0, waxman_smits_weight = 0, cec = -1, buckles_parameter = -1):
        """
        Calculates a petrophysical lithology and porosity model for conventional
        and unconventional reservoirs.

        For each depth, the method iterates a loop in 4 steps until convergence:

            1. Calculate Clay Volume

            Clay volume is calculated using a weighted method. Five different
            equations are available:

                I. Linear
                .. math::
                VCLAY = gr_index = \frac{GR_LOG - gr_matrix}{gr_clay - gr_matrix}

                II. Clavier
                .. math::
                VCLAY = 1.7 - \sqrt{3.38 - (gr_index + 0.7)^2}

                III. Larionov Tertiary Rocks
                .. math::
                VCLAY = 0.083 * (2^(3.7 * gr_index) - 1)

                IV. Neutron
                Calculate apparent nuetron log without organic matter, nphia.
                .. math::
                nphia = NPHI_LOG + (nphi_matrix - nphi_om) * vom

                Calculate vclay using neutron
                .. math::
                VCLAY = \frac{nphia - nphi_matrix}{nphi_clay - nphi_matrix}

                V. Neutron Density
                Calculate apprent density log without organic mater, rhoba.
                .. math::
                rhoba = RHOB_LOG + (rhom - rho_om) * vom

                Calculate vclay using neutron density
                .. math::
                m1 = \frac{nphi_fl - nphi_matrix}{rho_fl - rhom}
                x1 = nphi + m1 * (rhom - rhoba)
                x2 = nphi_clay + m1 * (rhom - rhoba)
                VCLAY = \frac{x1 - nphi_matrix}{x2 - nphi_matrix}

            First, the clay volume of the resepective euqations are calculated.
            They are then weighted with the vclay_weight. For example, if
            vclay_weight for every method is 1, then the final vclay is the
            average of the four equations. To use a single method, set
            vclay_method_weight = 1 and all other vclay_weight = 0.

            2. Calculate Total Organic Carbon And Pyrite

            TOC is calculated use a weighted method like vclay with three
            available equations:

                I. Schomker's Density Correlation
                .. math::
                TOC = schmoker_slope * (schmoker_baseline_rhob - RHOB_LOG)

                II. Passey's Nuetron Delta Log R
                .. math::
                dlr = 10^(\frac{RESDEEP_LOG}{passey_baseline_res} + 4 * (NPHI_LOG - passey_baseline_nphi)
                TOC = \frac{dlr_nphi * 10^(2.297 - 0.1688 * passey_lom)}{100}

                III. Passey's Density Delta Log R
                .. math::
                dlr = 10^(\frac{RESDEEP_LOG}{passey_baseline_res} - 2.5 * (RHOB_LOG - passey_baseline_rhob)
                TOC = \frac{dlr_nphi * 10^(2.297 - 0.1688 * passey_lom)}{100}

            For conventional reservoirs without organics, set vclay_cutoff = 1.

            3. Calculate Minerals And Porosity

            Non-negative least squares is used to find the remaining minerals
            according to the method described in Chapter 4 of Doveton's
            Principles of Mathematical Petrophysics.

            .. math::
            V = C^(-1) L

            Pure mineral log responses are required. For example, quartz would
            have the input:

                include_qtz = 'YES'
                rho_qtz = 2.65
                nphi_qtz = -0.04
                pe_qtz = 1.81

            An option to include exotic minerals is by specifiying the density,
            neutron, and pe response of mineral 'X'. For example, to add gypsum,
            use these parameters:

                include_x = 'YES'
                name_x = 'Gypsum'
                name_log_x = 'GYP'
                rho_x = 2.35
                nphi_x = 0.507
                pe_x = 4.04

            To exclude minerals because they are not present or essentially not
            present in the reservoir, set the include parameter to 'NO'. For
            example, to exclude dolomite:

                include_dol = 'NO'

            4. Calculate Saturations

            Saturation is calculated using a weighted method like vclay and toc.
            To use a single equation set equation_weight = 1 and all other
            equation_weight = 0. For example, to use only the Indonesia equation:

                archie_weight = 0
                simandoux_weight = 0
                modified_simandoux_weight = 0
                indonesia_weight = 1
                waxman_smits_weight = 0

            Five saturation equations are available:

                I. Archie
                .. math::
                SW = \frac{a * RW_LOG}{RESDEEP_LOG * phie^m}^(\frac{1}{n})

                II. Simandoux
                .. math::
                c = \frac{(1 - vclay) * a * RESDEEP_LOG}{phie^m}
                d = \frac{c * vclay}{2 * rt_clay}
                e = \frac{c}{RESDEEP_LOG}
                SW = ((d^2 + e)^2 - d)^(\frac{2}{n})

                III. Modified Simandoux
                .. math::

                IV. Indonesia (Poupon-Leveaux)
                .. math::
                f = \sqrt{\frac{phie^m}{RESDEEP_LOG}}
                g = \sqrt{\frac{vclay^(2 - vclay)}{rt_clay}}
                SW = ((f + g)^2 * RESDEEP_LOG)^(\frac{-1}{n})

                V. Waxman And Smits
                if cec <= 0:
                    .. math::
                    cec = 10^(1.9832 * vclay - 2.4473)
                else:
                    use input cec

                SW calculations:
                rw77 = RESDEEP_LOG * \frac{reservoir_temperature + 6.8}{83.8}
                b = 4.6 * (1 - 0.6 * e^(\frac{-0.77}{rw77}))
                f = a / phie^m
                qv = \frac{cec * (1 - phie) * rhom}{phie}
                SW = 0.5 * ((-b * qv * rw77) + \sqrt{(b * qv * rw77)^2 + \frac{4 * f * rw}{RESDEEP_LOG}})^(\frac{2}{n})


        Notes
        -----
        1. Clay bound water
            Clay bound water is included as part of the clay volume based on
            the default nphi_clay = 0.65. To calculate clay bound water
            seperately, set nphi_clay to clay matrix absent clay bound water,
            and include appropriate buckles_parameter for bound water
            saturations.

        2. Organics
            No differieniation is made between kerogen and other organic matter.


        Parameters
        ----------
        gr_matrix : float (default 10)
            Gamma Ray response of clean (non-clay) matrix
        nphi_matrix : float (default 0)
            Neutron response of clean (non-clay) matrix
        gr_clay : float (default 450)
            Gamma Ray response of pure clay matrix
        rho_clay : float (default 2.64)
            Density of pure clay matrix
        nphi_clay : float (default 0.65)
            Neutron reponse of pure clay matrix
        pe_clay : float (default 4)
            Photoelectric response of pure clay matrix
        rma : float (default 180)
            Resistivity of clean tight matrix
        rt_clay : float (default 80)
            Resistivity for inorganic shale matrix
        vclay_linear_weight : float (default 1)
            Weight of liner clay volume
        vclay_clavier_weight : float (defaul 0.5)
            Weight of Clavier clay volume
        vclay_larionov_weight : float (defaul 0.5)
            Weight of Larionov clay volume
        vclay_nphi_weight : float (default 1)
            Weight of Neutron clay volume
        vclay_nphi_rhob_weight : float (default 1)
            Weight of Neutron Density clay volume
        vclay_cutoff : float (default 0.05)
            Cutoff for oranics calculation. If vclay < vclay_cutoff then toc = 0
        rho_om : float (default 1.15)
            Density of organic matter
        nphi_om : float (default 0.6)
            Neutron response of pure organic matter
        pe_om : float (default 0.2)
            Photoelectric response of pure organic matter
        ro : float (default 1.6)
            Vitronite reflectance of organic matter
        lang_press : float (default 670)
            Langmiur pressure for gas adsorption on organic matter in psi
        passey_nphi_weight : float (default 1)
            Weight for Passey nphi toc
        passey_rhob_weight : float (default 1)
            Weight for Passey rhob toc
        passey_lom : float (default 10)
            Passey level of organic maturity
        passey_baseline_res : float (default 40)
            Passey inorganic baseline resistivity
        passey_baseline_rhob : float (default 2.65)
            Passey inorganic baseline density
        passey_baseline_nphi : float (default 0)
            Passey inorganic baseline neutron
        schmoker_weight : float (default 1)
            Weight for Schmoker toc
        schmoker_slope : float (default 0.7257)
            Slope for schmoker density to toc correlation
        schmoker_baseline_rhob : float (default 2.6)
            Density cutoff for schmoker toc correlation
        rho_pyr : float (default 5)
            Density of pyrite
        nphi_pyr : float (default 0.13)
            Neutron response of pure pyrite
        pe_pyr : float (default 13)
            Photoelectric response of pure pyrite
        om_pyrite_slope : float (default 0.2)
            Slope correlating pyrite volume to organic matter
        include_qtz : str {'YES', 'NO'} (default 'YES')
            Toggle to include or exclude qtz. 'YES' to include. 'NO' to exclude.
        rho_qtz : float (default 2.65)
            Density of quartz
        nphi_qtz : float (default -0.04)
            Neutron response for pure quartz
        pe_qtz : float (default 1.81)
            Photoelectric response for pure quartz
        include_clc : str {'YES', 'NO'} (default 'YES')
            Toggle to include or exclude clc. 'YES' to include. 'NO' to exclude.
        rho_clc : float (default 2.71)
            Density of calcite
        nphi_clc : float (default 0)
            Neutron response for pure calcite
        pe_clc : float (default 5.08)
            Photoelectric response for pure calcite
        include_dol : str {'YES', 'NO'} (default 'YES')
            Toggle to include or exclude dol. 'YES' to include. 'NO' to exclude.
        rho_dol : float (default 2.85)
            Density of dolomite
        nphi_dol : float (default 0.04)
            Neutron response to dolomite
        pe_dol : float (default 3.14)
            Photoelectric response to dolomite
        include_x : str {'YES', 'NO'} (default 'NO')
            Toggle to include or exclude exotic mineral, x. 'YES' to include.
            'NO' to exclude.
        name_x : str (default 'Gypsum')
            Name of exotic mineral, x.
        name_log_x : str (default 'GYP')
            Log name of exotic mineral, x
        rho_x : float (default 2.35)
            Density of exotic mineral, x
        nphi_x : float (default 0.507)
            Neutron response of exotic mineral, x
        pe_x : float (default 4.04)
            Photoelectric respone of exotic mineral, x
        pe_fl : float (default 0)
            Photoelectric response of reservoir fluid.
        m : float (default 2)
            Cementation exponent
        n : float (default 2)
            Saturation exponent
        a : float (default 1)
            Cementation constant
        cec : float (default -1)
            Cation Exchange Capaticy for use in Waxman Smits Sw equation. If
            cec = -1, correlation equation is used to calculate cec.
        archie_weight : float (default 0)
            Weight for archie Sw
        indonesia_weight : float (default 1)
            Weight for Indonesia Sw
        simandoux_weight : float (default 0)
            Weight for Simandoux Sw
        modified_simandoux_weight : float (default 0)
            Weight for Modified Simandoux Sw
        waxman_smits_weight : float (default 0)
            Weight for Waxman Smits Sw
        buckles_parameter : float (default -1)
            Buckles parameter for calculating irreducible water saturation. If
            less than 0, it is calculated using a correlation.

        Raises
        ------
        ValueError
            If fluid properties curve values are not present in self.curve_df,
            then ValueError is raised with incorrect curve requirements

        ValueError
            If raw curves GR_N, NPHI_N, RHOB_N, and RESDEEP_N are not present,
            then ValueError is raised with incorrect curve requirements as
            raw curve is either not present or precondtioning has not been run.

        ValueError
            If no formation value factor is found, then ValueError is raised to
            satisfy the calculation requirements.



        References
        ----------

        VCLAY
        Clavier, C., W. Hoyle, and D. Meunier, 1971a, Quantitative
            interpretation of thermal neutron decay time logs: Part I.
            Fundamentals and techniques: Journal of Petroleum Technology, 23,
            743–755
        Clavier, C., W. Hoyle, and D. Meunier, 1971b, Quantitative
            interpretation of thermal neutron decay time logs: Part II.
            Interpretation example, interpretation accuracy, and timelapse
            technique: Journal of Petroleum Technology, 23, 756–763.
        Larionov VV (1969). Borehole Radiometry: Moscow, U.S.S.R., Nedra.

        NEED NEUTRON, NEUTRON - DENSITY

        TOC
        Passey, Q. R., Creaney, S., Kulla, J. B., Moretti, F. J., Stroud, J. D.,
            1990, Practical Model for Organic Richness from Porosity and
            Resistivity Logs, AAPG Bulletin, 74, 1777-1794
        Schmoker, J.W., 1979, Determination of organic content of Appalachian
            Devonian shales from formation-density logs: American Association of
            Petroleum Geologists Bulletin, v.63, no.9, p.1504-1509

        MATRIX
        Doveton, John H. Principles of Mathematical Petrophysics. Oxford:
            Oxford University Press, 2014.

        SATURATIONS
        Archie, G.E. 1942. The Electrical Resistivity Log as an Aid in
            Determining Some Reservoir Characteristics. Trans. of AIME 146 (1):
            54-62.
        Bardon, C., and Pied, B.,1969, Formation water saturation in shaly
            sands: Society of Professional Well Log Analysts 10th Annual Logging
            Symposium Transactions: Paper Z,19 pp.
        Poupon, A. and Leveaux, J. 1971. Evaluation of Water Saturations in
            Shaly Formations. The Log Analyst 12 (4)
        Simandoux, P., 1963, Dielectricmeasurements on porous media application
            to the measurement of water saturations: study of the behaviour of
            argillaceous formations: Revue de l'Institut Francais du Petrole 18,
            Supplementary Issue, p. 193-215.
        Waxman, M.H., and L.J.M. Smits, Electrical Conductivity in Oil Bearing
            Shaly Sands, Society of Petroleum Engineers Journal, June,
            p.107-122, 1968.

        """

        ### check for requirements ###

        df_columns = self.curve_df.columns.values

        required_raw_curves = ['GR_N', 'NPHI_N', 'RHOB_N', 'RESDEEP_N']
        for curve in required_raw_curves:
            if curve not in df_columns:
                raise ValueError('Raw curve %s not found and is required for multimineral_model.' % curve)

        required_curves_from_fluid_properties = ['DEPTH', 'RW', 'RHO_HC',
                  'RHO_W', 'NPHI_HC', 'NPHI_W', 'RES_TEMP', 'NES', 'PORE_PRESS']
        for curve in required_curves_from_fluid_properties:
            if curve not in df_columns:
                raise ValueError('Fluid Properties curve %s not found. Run fluid_properties before multimineral_model.' % curve)

        all_required_curves = required_raw_curves + required_curves_from_fluid_properties

        if 'BO' not in df_columns and 'BG' not in df_columns:
            raise ValueError('Formation Volume Factor required for multimineral_model. Run fluid_properties first.')

        if 'BO' in df_columns:
            hc_class = 'OIL'
        else:
            hc_class = 'GAS'

        if 'PE_N' in df_columns:
            use_pe = True
        else:
            use_ps = False

        df = self.curve_df.copy()
        if top is not None:
            df = df[df.DEPTH >= top]
        else:
            top = self.curve_df.DEPTH.min()

        if bottom is not None:
            df = df[df.DEPTH <= bottom]
        else:
            bottom = self.curve_df.DEPTH.max()

        ### initialize minerals ###
        if include_qtz.upper()[0] == 'Y':
            include_qtz = True
        else:
            include_qtz = False

        if include_clc.upper()[0] == 'Y':
            include_clc = True
        else:
            include_clc = False

        if include_dol.upper()[0] == 'Y':
            include_dol = True
        else:
            include_dol = False

        if include_x.upper()[0] == 'Y':
            include_x = True
            name_log_x = name_log_x.upper()
        else:
            include_x = False

        ### calculations over depths ###
        df['sample_rate'] = df.DEPTH.diff()
        df.loc[0, 'sample_rate'] = df.iloc[1].DEPTH
        for i, row in df.iterrows():

            ### check for null values in data, skip if true ###
            if True in list(row[all_required_curves].isnull()): continue

            ### initial parameters to start iterations ###
            phie = 0.1
            rhom = 2.68
            rho_fl = 1
            nphi_fl = 1
            vom = 0

            bvqtz_prev = 1
            bvclc_prev = 1
            bvdol_prev = 1
            bvx_prev = 1
            phi_prev = 1
            bvom_prev = 1
            bvclay_prev = 1
            bvpyr_prev = 1

            diff = 1
            counter = 0
            while diff > 1 * 10 ** -3 and counter < 20:
                counter += 1

                ### log curves without organics ###
                rhoba = row.RHOB_N + (rhom - rho_om) * vom
                nphia = row.NPHI_N + (nphi_matrix - nphi_om) * vom

                ### clay solver ###
                gr_index = np.clip((row.GR_N - gr_matrix) / (gr_clay - gr_matrix), 0, 1)

                ### linear vclay method ###
                vclay_linear = gr_index

                ### Clavier vclay method ###
                vclay_clavier = np.clip(1.7 - np.sqrt(3.38 - (gr_index + 0.7) ** 2), 0, 1)

                ### larionov vclay method ###
                vclay_larionov = np.clip(0.083 * (2 ** (3.7 * gr_index) - 1), 0, 1)

                # Neutron vclay method without organic correction
                vclay_nphi = np.clip((nphia - nphi_matrix) / (nphi_clay - nphi_matrix), 0, 1)

                # Neutron Density vclay method with organic correction
                m1 = (nphi_fl - nphi_matrix) / (rho_fl - rhom)
                x1 = nphia + m1 * (rhom - rhoba)
                x2 = nphi_clay + m1 * (rhom - rho_clay)
                if x2 - nphi_matrix != 0:
                    vclay_nphi_rhob = np.clip((x1 - nphi_matrix) / (x2 - nphi_matrix), 0, 1)
                else:
                    vclay_nphi_rhob = 0

                vclay_weights_sum = vclay_linear_weight + vclay_clavier_weight + vclay_larionov_weight + vclay_nphi_weight + vclay_nphi_rhob_weight

                vclay = (vclay_linear_weight * vclay_linear + vclay_clavier_weight * vclay_clavier + vclay_larionov_weight * vclay_larionov + vclay_nphi_weight * vclay_nphi + vclay_nphi_rhob_weight * vclay_nphi_rhob) / vclay_weights_sum

                vclay = np.clip(vclay, 0, 1)

                bvclay = vclay * (1 - phie)

                ### organics ###
                if vclay > vclay_cutoff:

                    ### Passey, A Practical Model for Organic Richness from Porosity and Resistivity Logs, AAPG, 1990 Delta Log R ###
                    dlr_nphi = np.log10(row.RESDEEP_N / passey_baseline_res) + 4 * (row.NPHI_N - passey_baseline_nphi)
                    dlr_rhob = np.log10(row.RESDEEP_N / passey_baseline_res) - 2.5 * (row.RHOB_N - passey_baseline_rhob)

                    toc_nphi = np.clip((dlr_nphi * 10 ** (2.297 - 0.1688 * passey_lom) / 100), 0, 1)
                    toc_rhob = np.clip((dlr_rhob * 10 ** (2.297 - 0.1688 * passey_lom) / 100), 0, 1)

                    ### Schmoker ###
                    toc_sch = np.clip(schmoker_slope * (schmoker_baseline_rhob - row.RHOB_N), 0, 1)

                    toc_weights = passey_nphi_weight + passey_rhob_weight + schmoker_weight

                    ### toc in weight percent ###
                    toc = (passey_nphi_weight * toc_nphi + passey_rhob_weight * toc_rhob + schmoker_weight * toc_sch) / toc_weights

                    ### weight percent to volume percent ###
                    volume_om = toc / rho_om

                    rhom_no_om = (rhom - toc * rho_om) / (1 - toc) # matrix density without organic matter
                    volume_else = (1 - toc) / rhom_no_om

                    volume_total = volume_om + volume_else

                    vom = volume_om / volume_total
                    bvom = vom * (1 - phie)

                else:
                    toc = 0
                    vom = 0
                    bvom = 0

                ### pyrite correlation with organics ###
                vpyr = np.clip(om_pyrite_slope * vom, 0, 1)
                bvpyr = vpyr * (1 - phie)

                ### create C, V, and L matrix for equations in Chapter 4 of ###
                ### Principles of Mathematical Petrophysics by Doveton ###

                ### removed effect of clay, organics, and pyrite ###
                volume_unconventional = bvom + bvclay + bvpyr
                rhob_clean = (row.RHOB_N - (rho_om * bvom + rho_clay * bvclay + rho_pyr * bvpyr)) / (1 - volume_unconventional)
                nphi_clean = (row.NPHI_N - (nphi_om * bvom + nphi_clay * bvclay + nphi_pyr * bvpyr)) / (1 - volume_unconventional)

                minerals = []
                if use_pe:
                    pe_clean = (row.PE_N - (pe_om * bvom + pe_clay * bvclay + pe_pyr * bvpyr)) / (1 - bvom - bvclay - bvpyr)
                    l_clean = np.asarray([rhob_clean, nphi_clean, pe_clean, 1])
                    l = np.asarray([row.RHOB_N, row.NPHI_N, row.PE_N, 1])

                    c_clean = np.asarray([0,0,0]) # initialize mineral matrix C

                    if include_qtz:
                        minerals.append('QTZ')
                        mineral_matrix = np.asarray((rho_qtz, nphi_qtz, pe_qtz))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_clc:
                        minerals.append('CLC')
                        mineral_matrix = np.asarray((rho_clc, nphi_clc, pe_clc))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_dol:
                        minerals.append('DOL')
                        mineral_matrix = np.asarray((rho_dol, nphi_dol, pe_dol))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_x:
                        minerals.append('X')
                        mineral_matrix = np.asarray((rho_x, nphi_x, pe_x))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    fluid_matrix = np.asarray((rho_fl, nphi_fl, pe_fl))
                    c_clean = np.vstack((c_clean, fluid_matrix))
                    minerals.append('PHI')

                else:
                    l_clean = np.asarray([rhob_clean, nphi_clean, 1])
                    l = np.asarray([row.RHOB_N, row.NPHI_N, 1])

                    c_clean = np.asarray((0,0)) # initialize mineral matrix C

                    if include_qtz:
                        minerals.append('QTZ')
                        mineral_matrix = np.asarray((rho_qtz, nphi_qtz))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_clc:
                        minerals.append('CLC')
                        mineral_matrix = np.asarray((rho_clc, nphi_clc))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_dol:
                        minerals.append('DOL')
                        mineral_matrix = np.asarray((rho_dol, nphi_dol))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    if include_x:
                        minerals.append('X')
                        mineral_matrix = np.asarray((rho_x, nphi_x))
                        c_clean = np.vstack((c_clean, mineral_matrix))

                    fluid_matrix = np.asarray((rho_fl, nphi_fl))
                    c_clean = np.vstack((c_clean, fluid_matrix))
                    minerals.append('PHI')

                c_clean = np.delete(c_clean, 0, 0)

                c_clean = np.vstack((c_clean.T, np.ones_like(c_clean.T[0])))

                bv_clean = nnls(c_clean, l_clean.T)[0]

                bvqtz = 0
                bvclc = 0
                bvdol = 0
                bvx = 0

                component_sum = np.sum(bv_clean)

                for s, mineral in enumerate(minerals):
                    if mineral == 'QTZ':
                        bvqtz = (bv_clean[s] / component_sum) * (1 - volume_unconventional)
                        bv_clean[s] = bvqtz
                    if mineral == 'CLC':
                        bvclc = (bv_clean[s] / component_sum) * (1 - volume_unconventional)
                        bv_clean[s] = bvclc
                    if mineral == 'DOL':
                        bvdol = (bv_clean[s] / component_sum) * (1 - volume_unconventional)
                        bv_clean[s] = bvdol
                    if mineral == 'X':
                        bvx = (bv_clean[s] / component_sum) * (1 - volume_unconventional)
                        bv_clean[s] = bvx
                    if mineral == 'PHI':
                        phie = (bv_clean[s] / component_sum) * (1 - volume_unconventional)
                        bv_clean[s] = phie

                if use_pe:
                    c = np.hstack((c_clean, np.asarray( ((rho_om, rho_clay, rho_pyr), (nphi_om, nphi_clay, nphi_pyr), (pe_om, pe_clay, pe_pyr), (1, 1, 1)) )))
                else:
                    c = np.vstack((c_clean, np.asarray( ((rho_om, rho_clay, rho_pyr), (nphi_om, nphi_clay, nphi_pyr), (1, 1, 1)) )))

                bv = np.append(bv_clean, (bvom, bvclay, bvpyr))

                l_hat = np.dot(c, bv)

                sse = np.dot((l - l_hat).T, l - l_hat)

                prev = np.asarray((bvqtz_prev, bvclc_prev, bvdol_prev, bvx_prev, phi_prev, bvom_prev, bvclay_prev, bvpyr_prev))
                cur = np.asarray((bvqtz, bvclc, bvdol, bvx, phie, bvom, bvclay, bvpyr))

                diff = np.abs(cur - prev).sum()

                bvqtz_prev = bvqtz
                bvclc_prev = bvclc
                bvdol_prev = bvdol
                bvx_prev = bvx
                bvom_prev = bvom
                bvclay_prev = bvclay
                bvpyr_prev = bvpyr
                phi_prev = phie

                avg_percent_error = np.mean(np.abs(l - l_hat) / l) * 100

                ### calculate matrix volume fraction ###

                per_matrix = 1 - phie

                vqtz = bvqtz / per_matrix
                vclc = bvclc / per_matrix
                vdol = bvdol / per_matrix
                vx = bvx / per_matrix
                vclay = bvclay / per_matrix
                vom = bvom / per_matrix
                vpyr = bvpyr / per_matrix

        		### calculate weight fraction ###

                mass_qtz = vqtz * rho_qtz
                mass_clc = vclc * rho_clc
                mass_dol = vdol * rho_dol
                mass_x = vx * rho_x
                mass_om = vom * rho_om
                mass_clay = vclay * rho_clay
                mass_pyr = vpyr * rho_pyr

                rhom = mass_qtz + mass_clc + mass_dol + mass_x + mass_om + mass_clay + mass_pyr

                wtqtz = mass_qtz / rhom
                wtclc = mass_clc / rhom
                wtdol = mass_dol / rhom
                wtx = mass_x / rhom
                wtom = mass_om / rhom
                wtclay = mass_clay / rhom
                wtpyr = mass_pyr / rhom
                toc = wtom

                ### saturations ###

                ### porosity cutoff in case phie =  0 ###
                if phie < 0.001:
                    phis = 0.001
                else:
                    phis = phie

                ### Archie ###
                sw_archie = np.clip(((a * row.RW) / (row.RESDEEP_N * (phis ** m))) ** (1 / n), 0, 1)

                ### Indonesia ###
                sw_ind_a = (phie ** m / row.RW) ** 0.5
                sw_ind_b = (vclay ** (2.0 - vclay) / rt_clay) ** 0.5
                sw_indonesia = np.clip(((sw_ind_a + sw_ind_b) ** 2.0 * row.RESDEEP_N) ** (-1 / n), 0, 1)

                ### Simandoux ###
                c = (1.0 - vclay) * a * row.RW / (phis ** m)
                d = c * vclay / (2.0 * rt_clay)
                e = c / row.RESDEEP_N
                sw_simandoux = np.clip(((d**2 + e) ** 0.2 - d) ** (2 / n), 0, 1)

                ### modified Simandoux ###
                sw_modified_simandoux = np.clip((0.5 * row.RW / phis ** m) * ((4 * phis **m) / (row.RW * row.RESDEEP_N) + (vclay / rt_clay) ** 2) ** (1 / n) - vclay / rt_clay, 0, 1)

                ### Waxman Smits ###
                if cec <= 0:
                    cec = 10 ** (1.9832 * vclay - 2.4473)

                rw77 = row.RESDEEP_N * (row.RES_TEMP + 6.8) / 83.8
                b = 4.6 * (1 - 0.6 * np.exp(-0.77 / rw77))
                f = a / (phis ** m)
                qv = cec * (1 - phis) * rhom / phis
                sw_waxman_smits = np.clip(0.5 * ((-b * qv * rw77) + ((b * qv * rw77) ** 2 + 4 * f * row.RW / row.RESDEEP_N) ** 0.5 ) ** (2 / n), 0, 1)

                ### weighted calculation with bv output ###
                weight_saturations = archie_weight + indonesia_weight + simandoux_weight + modified_simandoux_weight + waxman_smits_weight

                sw = (archie_weight * sw_archie + indonesia_weight * sw_indonesia + simandoux_weight * sw_simandoux + modified_simandoux_weight * sw_modified_simandoux + waxman_smits_weight * sw_waxman_smits) / weight_saturations

                bvw = phie * sw
                bvh = phie * (1 - sw)

                if hc_class == 'OIL':
                    oip = (7758 * 640 * row.sample_rate * bvh * 10 ** -6) / row.BO # Mmbbl per sample rate

                elif hc_class == 'GAS':
                    langslope = (-0.08 * row.RES_TEMP + 2 * ro + 22.75) / 2
                    gas_ads = langslope * vom * 100 * (row.PORE_PRESS / (row.PORE_PRESS + lang_press))

                    gip_free = (43560 * 640 * row.sample_rate * bvh * 10 ** -9) / row.BG   # BCF per sample rate
                    gip_ads = (1359.7 * 640 * row.sample_rate * row.RHOB_N * gas_ads * 10 ** -9)/ row.BG	# BCF per sample rate
                    gip = gip_free + gip_ads

                rho_fl = row.RHO_W * sw + row.RHO_HC * (1 - sw)
                nphi_fl = row.NPHI_W * sw + row.NPHI_HC * (1 - sw)

            ### save calculations to dataframe ###

            ### bulk volume ###
            self.curve_df.loc[i, 'BVOM'] = bvom
            self.curve_df.loc[i, 'BVCLAY'] = bvclay
            self.curve_df.loc[i, 'BVPYR'] = bvpyr

            if include_qtz:
                self.curve_df.loc[i, 'BVQTZ'] = bvqtz
            if include_clc:
                self.curve_df.loc[i, 'BVCLC'] = bvclc
            if include_dol:
                self.curve_df.loc[i, 'BVDOL'] = bvdol
            if include_x:
                self.curve_df.loc[i,  'BV' + name_log_x] = bvx

            self.curve_df.loc[i, 'BVH'] = bvh
            self.curve_df.loc[i, 'BVW'] = bvw

            ### porosity and saturations ###
            self.curve_df.loc[i, 'PHIE'] = phie
            self.curve_df.loc[i, 'SW'] = sw
            self.curve_df.loc[i, 'SHC'] = 1 - sw

            ### mineral volumes ###
            self.curve_df.loc[i, 'VOM'] = vom
            self.curve_df.loc[i, 'VCLAY'] = vclay
            self.curve_df.loc[i, 'VPYR'] = vpyr

            if include_qtz:
                self.curve_df.loc[i, 'VQTZ'] = vqtz
            if include_clc:
                self.curve_df.loc[i, 'VCLC'] = vclc
            if include_dol:
                self.curve_df.loc[i, 'VDOL'] = vdol
            if include_x:
                self.curve_df.loc[i, 'V' + name_log_x] = vx

            ### weight percent ###
            self.curve_df.loc[i, 'RHOM'] = rhom
            self.curve_df.loc[i, 'TOC'] = toc
            self.curve_df.loc[i, 'WTCLAY'] = wtclay
            self.curve_df.loc[i, 'WTPYR'] = wtpyr

            if include_qtz:
                self.curve_df.loc[i, 'WTQTZ'] = wtqtz
            if include_clc:
                self.curve_df.loc[i, 'WTCLC'] = wtclc
            if include_dol:
                self.curve_df.loc[i, 'WTDOL'] = wtdol
            if include_x:
                self.curve_df.loc[i, 'WT' + name_log_x] = wtx

            ### find irreducible water if buckles_parameter is specified ###
            if buckles_parameter > 0:
                sw_irr = buckles_parameter / (phie / (1 - vclay))
                bvwi = phie * sw_irr
                bvwf = bvw - bvwi
                self.curve_df.loc[i, 'BVWI'] = bvwi
                self.curve_df.loc[i, 'BVWF'] = bvwf

            if hc_class == 'OIL':
                self.curve_df.loc[i, 'OIP'] = oip

            elif hc_class == 'GAS':
                self.curve_df.loc[i, 'GIP_FREE'] = gip_free
                self.curve_df.loc[i, 'GIP_ADS'] = gip_ads
                self.curve_df.loc[i, 'GIP'] = gip

        ### find irreducible water saturation outside of minerology loop since parameters depend on calculated values ###
        if buckles_parameter < 0:
            buckles_parameter = np.mean(self.curve_df[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)].PHIE * self.curve_df[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)].SW)

            sw_irr = buckles_parameter / (self.curve_df.PHIE / (1 - self.curve_df.VCLAY))
            self.curve_df.loc[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom), 'BVWI'] = self.curve_df[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)].PHIE * sw_irr
            self.curve_df.loc[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom), 'BVWF'] = self.curve_df[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)].BVW - self.curve_df[(self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)].BVWI

        ### append curves to curve_values and curve_parameters ###

        ### bulk volumes ###
        self.curve_parameters.append('BVOM')
        self.curve_values['BVOM'] = Parameter(name = 'BVOM', unit = 'v/v', des = 'Bulk Volume Fraction Organic Matter')

        self.curve_parameters.append('BVCLAY')
        self.curve_values['BVCLAY'] = Parameter(name = 'BVCLAY', unit = 'v/v', des = 'Bulk Volume Fraction Clay')

        self.curve_parameters.append('BVPYR')
        self.curve_values['BVPYR'] = Parameter(name = 'BVPYR', unit = 'v/v', des = 'Bulk Volume Fraction Pyrite')

        if include_qtz:
            self.curve_parameters.append('BVQTZ')
            self.curve_values['BVQTZ'] = Parameter(name = 'BVQTZ', unit = 'v/v', des = 'Bulk Volume Fraction Quartz')
        if include_clc:
            self.curve_parameters.append('BVCLC')
            self.curve_values['BVCLC'] = Parameter(name = 'BVCLC', unit = 'v/v', des = 'Bulk Volume Fraction Calcite')
        if include_dol:
            self.curve_parameters.append('BVDOL')
            self.curve_values['BVDOL'] = Parameter(name = 'BVDOL', unit = 'v/v', des = 'Bulk Volume Fraction Dolomite')
        if include_x:
            self.curve_parameters.append('BV' + name_log_x)
            self.curve_values['BV' + name_log_x] = Parameter(name = 'BV' + name_log_x, unit = 'v/v', des = 'Bulk Volume Fraction ' + name_x)

        self.curve_parameters.append('BVH')
        self.curve_values['BVH'] = Parameter(name = 'BVH', unit = 'v/v', des = 'Bulk Volume Fraction Hydrocarbon')

        self.curve_parameters.append('BVW')
        self.curve_values['BVW'] = Parameter(name = 'BVW', unit = 'v/v', des = 'Bulk Volume Fraction Water')

        self.curve_parameters.append('BVWI')
        self.curve_values['BVWI'] = Parameter(name = 'BVWI', unit = 'v/v', des = 'Bulk Volume Fraction Water irreducible')

        self.curve_parameters.append('BVWF')
        self.curve_values['BVWF'] = Parameter(name = 'BVWF', unit = 'v/v', des = 'Bulk Volume Fraction Water Free')

        ### porosity and saturations ###
        self.curve_parameters.append('PHIE')
        self.curve_values['PHIE'] = Parameter(name = 'PHIE', unit = 'v/v', des = 'Effective Porosity')

        self.curve_parameters.append('SW')
        self.curve_values['SW'] = Parameter(name = 'SW', unit = 'v/v', des = 'Water Saturation')

        self.curve_parameters.append('SHC')
        self.curve_values['SHC'] = Parameter(name = 'SHC', unit = 'v/v', des = 'Hydrocarbon Saturation Saturation')


        ### matrix mineral volumes ###
        self.curve_parameters.append('VOM')
        self.curve_values['VOM'] = Parameter(name = 'VOM', unit = 'v/v', des = 'Matrix Volume Fraction Organic Matter')

        self.curve_parameters.append('VCLAY')
        self.curve_values['VCLAY'] = Parameter(name = 'VCLAY', unit = 'v/v', des = 'Matrix Volume Fraction Clay')

        self.curve_parameters.append('VPYR')
        self.curve_values['VPYR'] = Parameter(name = 'VPYR', unit = 'v/v', des = 'Matrix Volume Fraction Pyrite')

        if include_qtz:
            self.curve_parameters.append('VQTZ')
            self.curve_values['VQTZ'] = Parameter(name = 'VQTZ', unit = 'v/v', des = 'Matrix Volume Fraction Quartz')
        if include_clc:
            self.curve_parameters.append('VCLC')
            self.curve_values['VCLC'] = Parameter(name = 'VCLC', unit = 'v/v', des = 'Matrix Volume Fraction Calcite')
        if include_dol:
            self.curve_parameters.append('VDOL')
            self.curve_values['VDOL'] = Parameter(name = 'VDOL', unit = 'v/v', des = 'Matrix Volume Fraction Dolomite')
        if include_x:
            self.curve_parameters.append('V' + name_log_x)
            self.curve_values['V' + name_log_x] = Parameter(name = 'V' + name_log_x, unit = 'v/v', des = 'Matrix Volume Fraction ' + name_x)


        ### matrix weight fraction ###
        self.curve_parameters.append('RHOM')
        self.curve_values['RHOM'] = Parameter(name = 'RHOM', unit = 'g/cc', des = 'Matrix Density')

        self.curve_parameters.append('TOC')
        self.curve_values['TOC'] = Parameter(name = 'TOC', unit = 'wt/wt', des = 'Matrix Weight Fraction Organic Matter')

        self.curve_parameters.append('WTCLAY')
        self.curve_values['WTCLAY'] = Parameter(name = 'WTCLAY', unit = 'wt/wt', des = 'Matrix Weight Fraction Clay')

        self.curve_parameters.append('WTPYR')
        self.curve_values['WTPYR'] = Parameter(name = 'WTPYR', unit = 'wt/wt', des = 'Matrix Weight Fraction Pyrite')

        if include_qtz:
            self.curve_parameters.append('WTQTZ')
            self.curve_values['WTQTZ'] = Parameter(name = 'WTQTZ', unit = 'wt/wt', des = 'Matrix Weight Fraction Quartz')

        if include_clc:
            self.curve_parameters.append('WTCLC')
            self.curve_values['WTCLC'] = Parameter(name = 'WTCLC', unit = 'wt/wt', des = 'Matrix Weight Fraction Calcite')

        if include_dol:
            self.curve_parameters.append('WTDOL')
            self.curve_values['WTDOL'] = Parameter(name = 'WTDOL', unit = 'wt/wt', des = 'Matrix Weight Fraction Dolomite')

        if include_x:
            self.curve_parameters.append('WT' + name_log_x)
            self.curve_values['WT' + name_log_x] = Parameter(name = 'WT' + name_log_x, unit = 'wt/wt', des = 'Matrix Weight Fraction ' + name_x)

        ### hydrocarbon in place ###
        if hc_class == 'OIL':
            self.curve_parameters.append('OIP')
            self.curve_values['OIP'] = Parameter(name = 'OIP', unit = 'Mmbbl / section', des = 'Oil in Place')

        elif hc_class == 'GAS':
            self.curve_parameters.append('GIP')
            self.curve_values['GIP'] = Parameter(name = 'GIP', unit = 'BCF / section', des = 'Gas in Place')

            self.curve_parameters.append('GIP_FREE')
            self.curve_values['GIP_FREE'] = Parameter(name = 'GIP_FREE', unit = 'BCF / section', des = 'Free Gas in Place')

            self.curve_parameters.append('GIP_ADS')
            self.curve_values['GIP_ADS'] = Parameter(name = 'GIP_ADS', unit = 'BCF / section', des = 'Adsorbed Gas in Place')


    def formation_multimineral_model(self, formations = [], parameter = 'default'):
        """
        Calculate multimineral model over formations with loaded paramters

        Parameters
        ----------
        formations : list (default [])
            list of formations to calculate fluid properties over
        parameter : str (default 'default')
            name of parameter to use for fluid properties parameter settings
            loaded in method fluid_properties_parameters_from_csv

        See Also
        --------
        multimineral_model
            calculates multimineral model from log data
        multimineral_parameters_from_csv
            loads multimineral model parameters
        tops_from_csv
            loads tops of log from csv

        """

        for form in formations:
            top = self.tops[form]
            bottom = self.next_formation_depth(form)

            params = self.multimineral_parameters[parameter]

            self.multimineral_model(top = top, bottom = bottom, **params)


    def add_pay_flag(self, formations, flag = 1, less_than_or_equal = [], greater_than_or_equal = [], name = ''):
        """
        Add Pay Flag based on curve cutoffs

        Parameters
        ----------
        formations : list
            list of formations, which must be found in preloaded tops
        flag : float
            Numeric value of pay flag. If interval meets pay requirments, then
            flag value. Else, 0.
        less_than_or_equal : list (of tuples (CURVE, value), default [])
            pay flag cutoff where interval is flaged if CURVE is less than or
            equal to value. Must be list of tuples
        greater_than_or_equal : list (of tuples (CURVE, value), default [])
            pay flag cutoff where interval is flaged if CURVE is greater than or
            equal to value. Must be list of tuples
        name : str (default '')
            End name of pay flag. Defaults to numeric PAY_FLAG_1 and increasing
            values as more pays flags are added.

        Example
        -------
        >>> # loads Wolfcamp and adds pay flag based on resistivity
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> gtoe = [('RESDEEP_N', 20)]
        >>> forms = ['WFMPA', 'WFMPB', 'WFMPC', 'WFMPD']
        >>> log.add_pay_flag(forms, greater_than_or_equal = gtoe, name = 'RES')

        """

        if len(name) < 1:
            c = 1
            for column in self.curve_df.columns.values:
                if 'PAY_FLAG' in column:
                    c += 1
            name = 'PAY_FLAG_' + str(c)
        else:
            name = 'PAY_FLAG_' + name

        for form in formations:
            top = self.tops[form]
            bottom = self.next_formation_depth(form)

            depths = (self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)
            cutoff = depths.copy()

            for curve, value in less_than_or_equal:
                cutoff = (cutoff) & (self.curve_df[depths][curve] <= value)

            for curve, value in greater_than_or_equal:
                cutoff = (cutoff) & (self.curve_df[depths][curve] >= value)

            self.curve_df.loc[depths, name] = np.nan
            self.curve_df.loc[cutoff, name] = flag
            self.curve_df.loc[(depths) & (self.curve_df[name].isnull()), name] = 0


    def summations(self, formations = [], curves = ['PHIE']):
        """
        Cumulative summations over formations for given curves.

        Parameters
        ---------
        formations : list (default [])
            list of formations to calculate summations
        curves : list (default ['PHIE'])
            list of curves to calculated cumulative summations. Values in list
            must be present in the curve_df dataframe.

        Example
        -------
        Sum Oil in Place for Wolfcamp A
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> log.tops_from_csv() # default tops
        >>> log.fluid_properties_parameters_from_csv() # load default parameters
        >>> log.formation_fluid_properties(formations = 'WFMPA') # run with defaults
        >>> log.multimineral_parameters_from_csv() # load default parameters
        >>> log.formation_multimineral_model(formation = 'WFMPA') # run with defaults
        >>> # run summations for oil in place over Wolfcamp A
        >>> log.summations(formations = ['WFMPA'], curves = ['OIP'])

        See Also
        --------
        statistics
            Calculates statistics over given formations and curves

        """

        hc_columns = ['OIP', 'GIP', 'GIP_ADS', 'GIP_FREE']

        sample_rate = self.curve_df.DEPTH.diff()
        sample_rate.loc[0] = self.curve_df.DEPTH.iloc[1]
        for f in formations:
            top = self.tops[f]
            bottom = self.next_formation_depth(f)

            depths = (self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)
            for c in curves:
                name = c + '_SUM'
                unit = self.curve_values[c].unit
                des = 'Summation of ' + c
                series = self.curve_df[depths][c]

                ### include sample rate in summation for non hydrocarbon columns ###
                if c not in hc_columns:
                    series = series * sample_rate[depths]

                self.curve_df.loc[depths, name] = series.ix[::-1].cumsum()[::-1]
                self.curve_parameters.append(name)
                self.curve_values[name] = Parameter(name = name, unit = unit, des = des)


    def statistics(self, formations = [], curves = ['PHIE']):
        """
        Curve statistcs for given formations and curves

        Paramters
        ---------
        formations : list (default [])
            list of formations to calculate statistics over. Must be included in
            preloaded tops
        curves : list (default ['PHIE'])
            list of curve to calculate statistics for. Must be includes in
            curves_df

        Returns
        -------
        df : DataFrame
            Returns Mean, Sum, and Standard Deviation for each curve over every
            formation

        Example
        -------
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> log.tops_from_csv() # default tops
        >>> f = ['WFMPA', 'WFMPB', 'WFMPC', 'WFMPD']
        >>> c = ['GR_N', 'RHOB_N', NPHI_N']
        >>> stats_df = log.statistics(formations = f, curves = c)
        >>> print(stats_df)

        COPY AND PASTE EXPORT DATA HERE!!!
        ###
        :LKSJDF:LSKJDFLKJS

        SEE ALSO
        --------
        statistics_to_csv
            saves statistics to csv file with option to append or overwrite
            formation values
        summations
            calculate summation curves over a given formation

        """

        hc_columns = ['OIP', 'GIP', 'GIP_ADS', 'GIP_FREE']

        sample_rate = self.curve_df.DEPTH.diff()
        sample_rate.loc[0] = self.curve_df.DEPTH.iloc[1]

        pay_flags = []
        for c in self.curve_df.columns.values:
            if 'PAY_FLAG' in c:
                pay_flags.append(c)

        stats_data = {}
        for f in formations:

            top = self.tops[f]
            bottom = self.next_formation_depth(f)
            formation_data = {'DATETIME': dt.datetime.now(), 'GROSS_H': bottom - top}

            depths = (self.curve_df.DEPTH >= top) & (self.curve_df.DEPTH < bottom)

            for c in curves:
                if c not in self.curve_df.columns.values:
                    raise ValueError('Column %s not in log curves.' % c)

                series = self.curve_df[depths][c]
                formation_data[c + '_MEAN'] = series.mean()
                formation_data[c + '_STD'] = series.std()

                if c not in hc_columns:
                    ### multiply by step rate for summations ###
                    series_sum = (series * sample_rate).sum()
                else:
                    series_sum = series.sum()

                formation_data[c + '_SUM'] = series_sum

                for p in pay_flags:
                    pay_series = series[(depths) & (self.curve_df[p] > 0)].copy()

                    formation_data[c + '_' + p + '_MEAN'] = pay_series.mean()
                    formation_data[c + '_' + p + '_STD'] = pay_series.std()

                    if c not in hc_columns:
                        ### multiply by step rate for summations ###
                        pay_series_sum = (pay_series * sample_rate).sum()
                    else:
                        pay_series_sum = pay_series.sum()

                    formation_data[c + '_' + p + '_SUM'] = pay_series_sum

            stats_data[f] = formation_data

        df = pd.DataFrame.from_dict(stats_data, orient = 'index')
        df.index.name = 'FORMATION'
        df['UWI'] = self.uwi
        df.reset_index(inplace = True)

        return df

    def statistics_to_csv(self, file_path, replace = False, formations = [], curves = ['PHIE']):
        """
        Saves curve statistcs for given formations and curves to a csv

        Paramters
        ---------
        file_path : str
            path to csv file
        replace : boolean (default False)
            option to replace uwi, formation statistics if already in csv file
        formations : list (default [])
            list of formations to calculate statistics over. Must be included in
            preloaded tops
        curves : list (default ['PHIE'])
            list of curve to calculate statistics for. Must be includes in
            curves_df

        Example
        -------
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> log.tops_from_csv() # default tops
        >>> f = ['WFMPA', 'WFMPB', 'WFMPC', 'WFMPD']
        >>> c = ['GR_N', 'RHOB_N', NPHI_N']
        >>> p = 'path/to/my/file.csv'
        >>> log.statistics_to_csv(p, formations = f, curves = c)

        SEE ALSO
        --------
        statistics
            calculates statistics and returns a dataframe
        summations
            calculate summation curves over a given formation

        """

        new_df = self.statistics(formations = formations, curves = curves)

        try:
            prev_df = pd.read_csv(file_path)
            prev_df['UWI'] = prev_df['UWI'].apply(str)
        except:
            prev_df = pd.DataFrame([])

        if replace and len(prev_df) > 0:
            for i, row in new_df.iterrows():
                drop_indexes = prev_df[(prev_df.UWI == row.UWI) & (prev_df.FORMATION == row.FORMATION)].index
                prev_df.drop(drop_indexes, inplace = True)
            new_df = prev_df.append(new_df)
        else:
            new_df = prev_df.append(new_df)

        new_df = new_df.set_index(['UWI', 'FORMATION'])

        new_df.to_csv(file_path)


    def write(self, file_path, top_depth = None, bottom_depth = None):
        """
        Write an las file of the log

        Parameters
        ----------
        file_path : str
            a path to save the las
        top_depth : float (default None)
            the start depth of the las file. If None, will begin at the top of
            recorded data
        bottom_depth : float (default None)
            end depth of the las file. If None, will export to the end of
            recorded data

        Example
        -------
        >>> # Read las file, then write to new location with better formatting
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> file_name = 'path/to/new/location/name_of_file.las'
        >>> log.write(file_name)

        >>> # Read las file, and only write Wolfcamp A for a smaller file size
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> log.tops_from_csv() # default tops
        >>> start = log.tops['WFMPA']
        >>> end = log.tops['WFMPB']
        >>> file_name = 'path/to/new/location/name_of_file.las'
        >>> log.write(file_name, top_depth = start, bottom_depth = end)

        """

        output = '~VERSION INFORMATION\n'
        for attribute in self.version_parameters:
            parameter = self.version_values[attribute]
            output += parameter.to_string() + '\n'

        output += '#'.ljust(64, '-') + '\n'
        output += '~WELL INFORMATION\n'
        output += '#MNEM'.ljust(12, ' ') + 'UNIT'.ljust(9, ' ') + 'DATA'.ljust(34, ' ') + 'DESCRIPTION\n'
        output += '#'.ljust(11, '-') + ' '.ljust(9, '-') + ' '.ljust(34, '-') + ' '.ljust(11 , '-') + '\n'
        for attribute in self.well_parameters:
            parameter = self.well_values[attribute]
            output += parameter.to_string() + '\n'

        output += '#'.ljust(64, '-') + '\n'
        output += '~PARAMETER INFORMATION\n'
        output += '#MNEM'.ljust(12, ' ') + 'UNIT'.ljust(8, ' ') + 'VALUE'.ljust(36, ' ') + 'DESCRIPTION\n'
        output += '#'.ljust(11, '-') + ' '.ljust(8, '-') + ' '.ljust(36, '-') + ' '.ljust(11 , '-') + '\n'
        for attribute in self.parameter_parameters:
            parameter = self.parameter_values[attribute]
            output += parameter.to_string() + '\n'

        output += '#'.ljust(64, '-') + '\n'
        output += '~CURVE INFORMATION\n'
        output += '#MNEM'.ljust(12, ' ') + 'UNIT'.ljust(8, ' ') + 'API CODE'.ljust(36, ' ') + 'DESCRIPTION\n'
        output += '#'.ljust(11, '-') + ' '.ljust(8, '-') + ' '.ljust(36, '-') + ' '.ljust(11 , '-') + '\n'
        df_output_columns = []
        for attribute in self.curve_parameters:
            if attribute not in df_output_columns:
                parameter = self.curve_values[attribute]
                df_output_columns.append(parameter.name)
                output += parameter.to_string() + '\n'

        output += '#'.ljust(64, '-') + '\n'
        output += '#\n'

        output += '#   '
        for name in df_output_columns:
            output += name.ljust(12, ' ')
        output += '\n'
        output += '#\n'
        output += '~A\n'

        if top_depth is None:
            top_depth = self.curve_df.DEPTH.min()
        if bottom_depth is None:
            bottom_depth = self.curve_df.DEPTH.max()

        self.curve_df.fillna(value = self.null, inplace = True)
        curve_data = self.curve_df[df_output_columns].to_records(index = False)
        curve_data = curve_data[(curve_data.DEPTH >= top_depth) & (curve_data.DEPTH <= bottom_depth)]

        columns = len(curve_data[0]) - 1
        column_lengths = {}
        for i in range(columns + 1):
            column_lengths[i] = 0

        for row in curve_data:
            output += '    '
            for j, value in enumerate(row):
                value = '%.4f' % value
                if j == columns:
                    output += value
                else:
                    value = value.ljust(12, ' ')
                    output += value
            output += '\n'

        with open(file_path, 'w') as f:
            f.write(output)


    def to_csv(self, **kwargs):
        """
        Write the log DataFrame to a comma-=separated values (csv) file.

        Calls pandas.DataFrame.to_csv which includes these parameters:

        Parameters
        ----------
        path_or_buf : str or file handle (default None)
            File path or object, if None is provided the result is returned as
            a string.
        sep : str (default ‘,’)
            Field delimiter for the output file.
        na_rep : str (default ‘’)
            Missing data representation
        float_format : str (default None)
            Format string for floating point numbers
        columns : sequence, optional
            Columns to write
        header : boolean or list of string (default True)
            Write out column names. If a list of string is given it is assumed
            to be aliases for the column names
        index : boolean (default True)
            Write row names (index)
        index_label : str or sequence, or False (default None)
            Column label for index column(s) if desired. If None is given, and
            header and index are True, then the index names are used. A sequence
            should be given if the DataFrame uses MultiIndex. If False do not
            print fields for index names. Use index_label=False for easier
            importing in R
        mode : str
            Python write mode, default ‘w’
        encoding : str, optional
            A string representing the encoding to use in the output file,
            defaults to ‘ascii’ on Python 2 and ‘utf-8’ on Python 3.
        compression : str, optional
            a string representing the compression to use in the output file,
            allowed values are ‘gzip’, ‘bz2’, ‘xz’, only used when the first
            argument is a filename
        line_terminator : str (default '\n')
            The newline character or character sequence to use in the output
            file
        quoting : constant from csv module, optional
            defaults to csv.QUOTE_MINIMAL. If you have set a float_format then
            floats are converted to strings and thus csv.QUOTE_NONNUMERIC will
            treat them as non-numeric
        quotechar : str with length 1 (default ‘”’)
            character used to quote fields
        doublequote : boolean (default True)
            Control quoting of quotechar inside a field
        escapechar : str with length 1 (default None)
            character used to escape sep and quotechar when appropriate
        chunksize : int (default None)
            rows to write at a time
        tupleize_cols : boolean (default False)
            write multi_index columns as a list of tuples (if True) or new
            (expanded format) if False)
        date_format : str (default None)
            Format string for datetime objects
        decimal: str (default ‘.’)
            Character recognized as decimal separator. E.g. use ‘,’ for
            European data

        Example
        -------
        >>> # Read las file, then write to csv for use in excel
        >>> import petropy as ptr
        >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
        >>> file_name = 'path/to/save/location/name_of_file.csv'
        >>> log.to_csv(path_or_buf = file_name, index = False)

        """

        self.curve_df.fillna(value = self.null, inplace = True)
        self.curve_df.to_csv(**kwargs)
