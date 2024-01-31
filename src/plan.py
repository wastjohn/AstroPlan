# This python notebook shows the position of the Sun and also of a user-specified
# object.  --- JMC, Jan 19-22, 2024

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

# %matplotlib notebook
def plan(ra, dec, name):

    # Enter source ra,dec (J2000) in DEGREES:
    source=SkyCoord(ra*u.deg,dec*u.deg)

    # Sun ephemeris:
    time = Time('2024-01-01 00:00:00')
    times = time + np.linspace(0, 364, 365)*u.day
    sun = get_sun(times)

    ra = coord.Angle(ra*u.degree)
    dec = coord.Angle(dec*u.degree)

    fig, ax = plt.subplots(figsize=(15, 11))
    ax.set_xlabel('Right Ascension [deg]',size=30)
    ax.set_ylabel('Declination [deg]',size=30)
    ax.set_title(name,size=30)

    ax.tick_params(axis='both', labelsize=25)
    ax.set_xlim(0.*u.deg, 360*u.deg)
    # ax.set_ylim(-30*u.deg, 90*u.deg)
    ax.scatter(sun.ra, sun.dec, c='grey', linewidth=0, marker='.')

    ax.scatter(ra, dec, color='black', linewidth=1, s=128, label=str(name),marker='*')

    # Now color code by month:    
    janstarttime = Time('2024-01-01 12:00:00')
    jantimes = janstarttime + np.linspace(0, 30, 31)*u.day
    jansun = get_sun(jantimes)
    ax.scatter(jansun.ra, jansun.dec, c='purple', linewidth=2, marker='.')
    janone = get_sun(janstarttime)
    ax.scatter(janone.ra, janone.dec, c='purple',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Jan 1 2024')

    febstarttime = Time('2024-02-01 12:00:00')
    febtimes = febstarttime + np.linspace(0, 27, 28)*u.day
    febsun = get_sun(febtimes)
    ax.scatter(febsun.ra, febsun.dec, c='magenta', linewidth=2, marker='.')
    febone = get_sun(febstarttime)
    ax.scatter(febone.ra, febone.dec, c='magenta',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Feb 1 2024')

    marstarttime = Time('2024-03-01 12:00:00')
    martimes = marstarttime + np.linspace(0, 30, 31)*u.day
    marsun = get_sun(martimes)
    ax.scatter(marsun.ra, marsun.dec, c='pink', linewidth=2, marker='.')
    marone = get_sun(marstarttime)
    ax.scatter(marone.ra, marone.dec, c='pink',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Mar 1 2024')

    aprstarttime = Time('2024-04-01 12:00:00')
    aprtimes = aprstarttime + np.linspace(0, 29, 30)*u.day
    aprsun = get_sun(aprtimes)
    ax.scatter(aprsun.ra, aprsun.dec, c='cyan', linewidth=2, marker='.')
    aprone = get_sun(aprstarttime)
    ax.scatter(aprone.ra, aprone.dec, c='cyan',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Apr 1 2024')

    maystarttime = Time('2024-05-01 12:00:00')
    maytimes = maystarttime + np.linspace(0, 30, 31)*u.day
    maysun = get_sun(maytimes)
    ax.scatter(maysun.ra, maysun.dec, c='brown', linewidth=2, marker='.')
    mayone = get_sun(maystarttime)
    ax.scatter(mayone.ra, mayone.dec, c='brown',s=600,linewidth=2, edgecolor='black', marker='.', label = 'May 1 2024')

    junstarttime = Time('2024-06-01 12:00:00')
    juntimes = junstarttime + np.linspace(0, 29, 30)*u.day
    junsun = get_sun(juntimes)
    ax.scatter(junsun.ra, junsun.dec, c='black', linewidth=2, marker='.')
    junone = get_sun(junstarttime)
    ax.scatter(junone.ra, junone.dec, c='black',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Jun 1 2024')

    julstarttime = Time('2024-07-01 12:00:00')
    jultimes = julstarttime + np.linspace(0, 30, 31)*u.day
    julsun = get_sun(jultimes)
    ax.scatter(julsun.ra, julsun.dec, c='grey', linewidth=2, marker='.')
    julone = get_sun(julstarttime)
    ax.scatter(julone.ra, julone.dec, c='grey',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Jul 1 2024')

    augstarttime = Time('2024-08-01 12:00:00')
    augtimes = augstarttime + np.linspace(0, 30, 31)*u.day
    augsun = get_sun(augtimes)
    ax.scatter(augsun.ra, augsun.dec, c='red', linewidth=2, marker='.')
    augone = get_sun(augstarttime)
    ax.scatter(augone.ra, augone.dec, c='red',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Aug 1 2024')

    sepstarttime = Time('2024-09-01 12:00:00')
    septimes = sepstarttime + np.linspace(0, 29, 30)*u.day
    sepsun = get_sun(septimes)
    ax.scatter(sepsun.ra, sepsun.dec, c='orange', linewidth=2, marker='.')
    sepone = get_sun(sepstarttime)
    ax.scatter(sepone.ra, sepone.dec, c='orange',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Sep 1 2024')

    octstarttime = Time('2024-10-01 12:00:00')
    octtimes = octstarttime + np.linspace(0, 30, 31)*u.day
    octsun = get_sun(octtimes)
    ax.scatter(octsun.ra, octsun.dec, c='yellow', linewidth=2, marker='.')
    octone = get_sun(octstarttime)
    ax.scatter(octone.ra, octone.dec, c='yellow',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Oct 1 2024')

    novstarttime = Time('2024-11-01 12:00:00')
    novtimes = novstarttime + np.linspace(0, 29, 30)*u.day
    novsun = get_sun(novtimes)
    ax.scatter(novsun.ra, novsun.dec, c='limegreen', linewidth=2, marker='.')
    novone = get_sun(novstarttime)
    ax.scatter(novone.ra, novone.dec, c='limegreen',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Nov 1 2024')

    decstarttime = Time('2024-12-01 12:00:00')
    dectimes = decstarttime + np.linspace(0, 30, 31)*u.day
    decsun = get_sun(dectimes)
    ax.scatter(decsun.ra, decsun.dec, c='blue', linewidth=2, marker='.')
    decone = get_sun(decstarttime)
    ax.scatter(decone.ra, decone.dec, c='blue',s=600,linewidth=2, edgecolor='black', marker='.', label = 'Dec 1 2024')


    ax.legend(loc='upper right', fontsize=14,facecolor='grey', framealpha=0.2, edgecolor='black')
    return fig
# fig.savefig(outdir+str(name)+'.Plan-Source-Sun.png',dpi=300, bbox_inches = "tight", transparent = False, facecolor = 'whitesmoke')


# # Macalester information:
def getMacInfo():
    mac = Observer(longitude=-93.1691*u.deg, latitude=44.9379*u.deg, elevation=240*u.m, name="Macalester", timezone="US/Central")
    now=Time.now()
    sun= get_body('sun',now)
    deltasun = (sun.transform_to(AltAz(obstime=now,location=mac.location)).alt).value - (sun.transform_to(AltAz(obstime=now+1*u.second,location=mac.location)).alt).value
    if deltasun > 0:
        sunsetting = "Sun is setting"
    else:
        sunsetting = "Sun is rising"
    
    output = f"""
    ---------------------------------------
    Macalester is located at:\n
    \tLatitude = {str(mac.latitude)}\n
    \tLongitude = {str(mac.longitude)}\n
    \tAltitude = {str(mac.elevation)}\n
    \n
    Current Macalester information (UTC):\n
    \tDate = {str(now)} UTC\n
    \tJD   = {str(now.jd)}\n
    \tMJD  = {str(now.mjd)}\n
    \tDecimal year = {str(now.decimalyear)}\n
    \tLST = {str(mac.local_sidereal_time(now))}\n
    \tSun altitude = {str(sun.transform_to(AltAz(obstime=now,location=mac.location)).alt)}\n
    \t{sunsetting}\n
    \n
    Sunset and sunrise times (UTC):\n
    \tSunset = {str(mac.sun_set_time(now, which='nearest').to_datetime(mac.timezone))}\n
    \tMidnight = {str(mac.midnight(now, which='next').to_datetime(mac.timezone))}\n
    \tSunrise = {str(mac.sun_rise_time(now, which='next').to_datetime(mac.timezone))}\n
    ---------------------------------------
    """
    return output

# # RLMT information:
def getRLMTInfo():
    winer = Observer.at_site('Winer')
    now=Time.now()
    sun= get_body('sun',now)
    deltasun = (sun.transform_to(AltAz(obstime=now,location=winer.location)).alt).value - (sun.transform_to(AltAz(obstime=now+1*u.second,location=winer.location)).alt).value
    if deltasun > 0:
        sunsetting = "Sun is setting"
    else:
        sunsetting = "Sun is rising"
    
    output = f"""
    ---------------------------------------
    Winer is located at:\n
    \tLatitude = {str(winer.latitude)}\n
    \tLongitude = {str(winer.longitude)}\n
    \tAltitude = {str(winer.elevation)}\n
    \n
    Current Winer information (UTC):\n
    \tDate = {str(now)} UTC\n
    \tJD   = {str(now.jd)}\n
    \tMJD  = {str(now.mjd)}\n
    \tDecimal year = {str(now.decimalyear)}\n
    \tLST = {str(winer.local_sidereal_time(now))}\n
    \tSun altitude = {str(sun.transform_to(AltAz(obstime=now,location=winer.location)).alt)}\n
    \t{sunsetting}\n
    \n
    Sunset and sunrise times (UTC):\n
    \tSunset = {str(winer.sun_set_time(now, which='nearest').to_datetime(winer.timezone))}\n
    \tMidnight = {str(winer.midnight(now, which='next').to_datetime(winer.timezone))}\n
    \tSunrise = {str(winer.sun_rise_time(now, which='next').to_datetime(winer.timezone))}\n
    ---------------------------------------
    """
    return output


# Macalester Observatory

# Set cadence of plots to examine location of source, Sun, and Moon:
def getMacObsPlots(ra, dec, name):
    source=SkyCoord(ra*u.deg,dec*u.deg)
    mac = Observer(longitude=-93.1691*u.deg, latitude=44.9379*u.deg, elevation=240*u.m, name="Macalester", timezone="US/Central")

    starttime = Time.now()
    obs=EarthLocation(mac.location)
    interval = 7*u.day
    # interval = 1*u.day

    midnightlist = []

    for i in range(0,8):  # number of weeks in a year + 1
    # for i in range(0,365):
        time = starttime + interval*[i]
        midnight = (mac.midnight(time, which='next')).to_value('iso',subfmt='date_hm')
        date = time.to_value('iso',subfmt='date_hm')[0]
        midnighttime = Time(midnight)
        midnightlist.append(midnighttime)

    figs = []
    # Loop over midnightlist to plot the Sun, Moon, and source altitude vs. time:
    # with PdfPages(outdir+str(name)+'.Mac-weekly-plan.pdf') as pdf:
    for i in range(0,len(midnightlist)):
        time = midnightlist[i]    
        sourcealtaz = source.transform_to(AltAz(obstime=time,location=obs))
        delta_t = np.linspace(-12, 12, 2001)*u.hour
        frame = AltAz(obstime=time+delta_t,location=obs)
        sourcealtaz = source.transform_to(frame)
        times = time + delta_t
        sunaltaz = get_body("sun",times).transform_to(frame)
        moonaltaz = get_body("moon",times).transform_to(frame)

        fig, ax = plt.subplots(figsize=(15, 11))
        ax.set_xlabel(r'Hours from astronomical midnight at UT = '+str(time))
        ax.set_ylabel('Altitude [deg]',size=30)
        ax.tick_params(axis='both', labelsize=25)
        ax.set_xlim(-12.*u.hour, 12*u.hour)
        ax.set_ylim(0*u.deg, 90*u.deg)
        ax.plot(delta_t, sunaltaz.alt, color='orange', label='Sun', linewidth=3)
        ax.plot(delta_t, moonaltaz.alt, color='darkgrey', ls='--', label='Moon', linewidth=2)   
        ax.plot(delta_t, sourcealtaz.alt, color='dodgerblue', label=str(name), linewidth=4)  
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -6*u.deg, color='0.5', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -12*u.deg, color='0.2', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -18*u.deg, color='k', zorder=0)
        plt.legend(loc=1,fontsize=20, facecolor='white', framealpha=0.95)
        plt.title('Macalester Observatory: '+str(name)+' on '+str(midnightlist[i])[0:10],fontsize=25)
        # pdf.savefig(dpi=100)
        plt.close()
        figs.append(fig)
    
    return figs

# RLMT:
# Set cadence of plots to examine location of source, Sun, and Moon:

def getRLMTObsPlots(ra, dec, name):
    source=SkyCoord(ra*u.deg,dec*u.deg)
    winer = Observer.at_site('Winer')


    starttime = Time.now()
    obs=EarthLocation(winer.location)
    interval = 7*u.day
    # interval = 1*u.day

    midnightlist = []

    for i in range(0,8):  # number of weeks in a year + 1
    # for i in range(0,365):
        time = starttime + interval*[i]
        midnight = (winer.midnight(time, which='next')).to_value('iso',subfmt='date_hm')
        date = time.to_value('iso',subfmt='date_hm')[0]
        midnighttime = Time(midnight)
        midnightlist.append(midnighttime)

    figs = []
    # Loop over midnightlist to plot the Sun, Moon, and source altitude vs. time:
    # with PdfPages(outdir+str(name)+'.RLMT-weekly-plan.pdf') as pdf:
    for i in range(0,len(midnightlist)):
        time = midnightlist[i]    
        sourcealtaz = source.transform_to(AltAz(obstime=time,location=obs))
        delta_t = np.linspace(-12, 12, 2001)*u.hour
        frame = AltAz(obstime=time+delta_t,location=obs)
        sourcealtaz = source.transform_to(frame)
        times = time + delta_t
        sunaltaz = get_body("sun",times).transform_to(frame)
        moonaltaz = get_body("moon",times).transform_to(frame)

        fig, ax = plt.subplots(figsize=(15, 11))
        ax.set_xlabel(r'Hours from astronomical midnight at UT = '+str(time))
        ax.set_ylabel('Altitude [deg]',size=30)
        ax.tick_params(axis='both', labelsize=25)
        ax.set_xlim(-12.*u.hour, 12*u.hour)
        ax.set_ylim(0*u.deg, 90*u.deg)
        ax.plot(delta_t, sunaltaz.alt, color='orange', label='Sun', linewidth=3)
        ax.plot(delta_t, moonaltaz.alt, color='darkgrey', ls='--', label='Moon', linewidth=2)   
        ax.plot(delta_t, sourcealtaz.alt, color='dodgerblue', label=str(name), linewidth=4)  
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -6*u.deg, color='0.5', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -12*u.deg, color='0.2', zorder=0)
        ax.fill_between(delta_t, 0*u.deg, 90*u.deg, sunaltaz.alt < -18*u.deg, color='k', zorder=0)
        plt.legend(loc=1,fontsize=20, facecolor='white', framealpha=0.95)
        plt.title('RLMT: '+str(name)+' on '+str(midnightlist[i])[0:10],fontsize=25)
        # pdf.savefig(dpi=100)
        plt.close()
        figs.append(fig)
    return figs

def generateRLMTSchedule(title, observer, epoch, target, ra_hms, dec_dms, positioning, filter, duration, cmosmod, utstart, cadence, repeat):
    result = f"""
    title '{title}'\n
    observer '{observer}'\n
    epoch '{epoch}'\n

    source '{target}' ra {ra_hms} dec {dec_dms} positioning {positioning}\n
    filter {filter} duration {duration} cmosmod {cmosmod} utstart {utstart} cadence {cadence} repeat {repeat} /\n
    """
    return result