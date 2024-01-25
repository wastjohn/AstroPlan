"""
File: main.py

The goal of this file is to create a web app with steamlit to help identify and plan observation sessions with an optical telescope given a location and date, or a location and object.

Author: Will St. John
Date: January 2024
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astroplan import Observer
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from datetime import datetime
from astropy.coordinates import solar_system_ephemeris, FK5
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from astropy.coordinates import get_body
from astropy.table import Table
import astropy.coordinates as coord
from astropy.coordinates import angular_separation
from matplotlib.backends.backend_pdf import PdfPages
import plan

# %matplotlib notebook
st.title("AstroPlan")

st.write("This is a web app to help you plan your astrophotography sessions.")

"""
Let's start with my rational: I am in an observational astronomy course that requires I complete a photometry project using two telescopes. I will be observing star clusters. Before I can make any observation, I need to select an object. The objects I can see are limited by: the weather, sun position, moon position, and location of the telescope. 

I want to make a telescope dashboard that will help me plan an observation session. Ideally there are two possible ways this site could be used:
1. I know the type of object I want to observe, and I want to know what is/will be availble given my location and the date/dates I select
2. I know the location and date I want to obseve, and I want to see what objects are available to me.

For either option, I will need to have a database of objects I can query with RA's and Dec's.

An object is visible if:
- The object is above the horizon
- The sun is below the horizon
- The moon is not too close to the object
- The weather is good

After an option is selected, a series of plots will be generated based on the user's inputted coordinates and date. The plots will be:
- Sky view from the location (looking up)
- Moon's transit over time

If the user selects an object type, the database will be filtered to objects of that type. The transit times for each of them will be arranged by longest duration to shortest.

I will be using the following packages:
- streamlit
- astroplan

Website Proceedure:
- Get location, date, and object (if applicable, only need location --> use current night as default)

- If only location and date:
    - Make sun plot
    - Make moon plot
    - Make weather plot

- If location and object:
    - Make sun plot
    - make moon plot
    - make weather plot
    - make object plot



Ok. New goal: Turn John's code into a web app.
"""
st.write(plan.plan())