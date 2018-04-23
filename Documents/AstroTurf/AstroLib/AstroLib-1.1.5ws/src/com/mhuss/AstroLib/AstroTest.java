/** ***************************************************************************\
 * AstroTest
 *
 * A Generic Test Harness
 *
 * \**************************************************************************** */
package com.mhuss.AstroLib;

//import java.util.Date;
//import java.util.Calendar;
//import java.util.GregorianCalendar;
public class AstroTest
{
    public static void main(String args[])
    {
        // Punit();
        LunarRS();
    }

    public static void LunarRS()
    {
        // Long E is positive
        ObsInfo oi = new ObsInfo(new Latitude(42.4304), new Longitude(19.2594));
        System.out.println(LunarCalc.summary(oi));
        System.out.println(LunarCalc.summaryPHL());
    }

    public static void Punit()
    {
        double jd = DateOps.dmyToDay(11, 4, 1979);
        ObsInfo oi = new ObsInfo(new Latitude(27.20), new Longitude(77.02));

        // Sun
        PlanetData pde = new PlanetData(Planets.SUN, jd, oi);
        try
        {
            System.out.println("Sun Lon = " + Math.toDegrees(pde.getEclipticLon()));
        }
        catch (NoInitException e)
        {
            System.out.println("Error calculating Sun: " + e);
        }

        // Mars
        PlanetData pdm = new PlanetData(Planets.MARS, jd, oi);
        try
        {
            System.out.println("Mars Lon = " + Math.toDegrees(pdm.getEclipticLon()));
        }
        catch (NoInitException e)
        {
            System.out.println("Error calculating Mars: " + e);
        }
    }
}
