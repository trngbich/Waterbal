# -*- coding: utf-8 -*-
"""
Authors: Elga Salvadore
         IHE Delft 2019
Contact: e.salvadore@un-ihe.org
Repository: 
Module: WAWB


Description:
Simple water balance tool that uses RS data to compute major water balcen components at pixel level.
The code is a simplified version of WaterPix model (Authors: Gonzalo Espinoza and Claire Michailovsky)

"""
from .main import run

__all__ = ['run']

__version__ = '0.1'
