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
ra_col, dec_col, target_name_col = st.columns(3)
ra = ra_col.number_input("RA", min_value=0.0000000, max_value=359.9999999, value=83.8186621, step=0.0000001, format="%.7f")
dec = dec_col.number_input("DEC", min_value=-90.0000000, max_value=90.0000000, value=-5.3896789, step=0.0000001, format="%.7f")
target_name = target_name_col.text_input("Target Name", value="Orion Nebula")


st.write("Object location against Sun Emphemris")
st.write(plan.plan(ra,dec,target_name))
"""
How to read this chart:\n
The curve corresponds to the Sun's apparent position in the sky for a given time of the year.\n
The black star corresponds to the object's position in the sky.\n
If RA of the object is near the RA of the Sun for the given year, you won't be able to observe the object (because the object will be up when the Sun is up)!\n
Depending on your location (longitude and latitude), the object might be above or below the horizon.\n
Your latitude is the Dec of your zenith (the angle normal to the ground at your location)
So if your latitude is 40 degrees, your zenith is at 40 degrees Dec.\n
Assuming you can see an entire hemisphere of the sky, your horizon is at 40-90=-50 degrees and 40+90=130 degrees (max is 90), so you could objects with Decs ranging between -50 and 90 Dec.\n
"""
mac_col, rlmt_col = st.columns(2)

mac_col.write(plan.getMacInfo())
rlmt_col.write(plan.getRLMTInfo())

for im in plan.getMacObsPlots(ra,dec,target_name):
    mac_col.write(im)

for im in plan.getRLMTObsPlots(ra,dec,target_name):
    rlmt_col.write(im)
