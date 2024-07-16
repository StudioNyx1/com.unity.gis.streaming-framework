using System;
using System.Runtime.CompilerServices;

using Unity.Burst;
using UnityEngine;
using Unity.Mathematics;
using Unity.Geospatial.HighPrecision;

namespace Unity.Geospatial.Streaming
{
    /// <summary>
    /// Methods allowing to convert from / to <see cref="GeodeticCoordinates"/> when based on the Earth-fixed terrestrial
    /// reference system and geodetic datum: WGS84.
    /// </summary>
    public static class Wgs84
    {
        private const double k_a = 6_378_137.0;

        private const double k_b = (1 - k_f) * k_a;

        private const double k_f = 1 / 298.257223563;

        private const double k_ff = (1.0 - k_f) * (1.0 - k_f);

        private const double k_e2 = 1 - (k_b * k_b)/(k_a * k_a);
        
        private const double k_HALFE2 = k_e2 / 2;
        private const double k_P1ME2 = (1 - k_e2);
        private const double k_A2 = k_a * k_a;
        private const double k_INVA2 = 1.0 / k_A2;

        internal static double MajorRadius
        {
            get { return k_a; }
        }

        private static double MinorRadius
        {
            get { return k_b; }
        }

        private static double3 GeodeticToXzyEcef(double latitude, double longitude, double elevation)
        {
            return GeodeticToXzyEcef(new GeodeticCoordinates(latitude, longitude, elevation));
        }

        internal static double3 GeodeticToXzyEcef(GeodeticCoordinates coords)
        {
            TrigonometricRatios ratios = new TrigonometricRatios(coords);
            return GeodeticToXzyEcef(in ratios, coords.Elevation);
        }

        [BurstCompile(CompileSynchronously = true)]
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double3 GeodeticToXzyEcef(in TrigonometricRatios ratios, double elevation)
        {
            double c = 1.0 / math.sqrt(Square(ratios.CosLat) + k_ff * Square(ratios.SinLat));

            double s = c * k_ff;

            return new double3(
                (k_a * c + elevation) * ratios.CosLat * ratios.CosLon,
                (k_a * s + elevation) * ratios.SinLat,
                (k_a * c + elevation) * ratios.CosLat * ratios.SinLon);
        }

        /// <summary>
        /// Convert a geodetic (latitude/longitude) to a euclidean transformation.
        /// </summary>
        /// <param name="position">The coordinates values representing the translation expressed in degrees and minutes..</param>
        /// <param name="eulerAngles">Orientation of the vector where zero (0, 0, 0) is pointing to the same direction as the <see cref="Position"/> normal.</param>
        /// <returns>The converted translation / rotation.</returns>
        public static EuclideanTR GeodeticToXzyEcef(GeodeticCoordinates position, float3 eulerAngles)
        {
            EuclideanTR result;

            double4x4 newMatrix = GetXzyEcefFromXzyEnuMatrix(position);

            newMatrix.GetTRS(out result.Position, out quaternion geodeticIdentityRotation, out _);

            result.Rotation = math.mul(geodeticIdentityRotation, Quaternion.Euler(FlipPrincipalAxes(eulerAngles)));

            return result;
        }

        /// <summary>
        /// Convert a geodetic (latitude/longitude) to a euclidean transformation.
        /// </summary>
        /// <param name="position">The coordinates values representing the translation expressed in degrees and minutes..</param>
        /// <param name="rotation">The rotation where zero (0, 0, 0) is pointing to the same direction as the <see cref="Position"/> normal.</param>
        /// <returns>The converted translation / rotation.</returns>
        public static EuclideanTR GeodeticToXzyEcef(GeodeticCoordinates position, quaternion rotation)
        {
            EuclideanTR result;

            double4x4 newMatrix = GetXzyEcefFromXzyEnuMatrix(position);

            newMatrix.GetTRS(out result.Position, out quaternion geodeticIdentityRotation, out _);

            FlipPrincipalAxes(ref rotation);

            result.Rotation = math.mul(geodeticIdentityRotation, rotation);

            return result;
        }

        /// <summary>
        /// Convert a geodetic (latitude/longitude) to a euclidean matrix.
        /// </summary>
        /// <param name="origin">The position to convert as a matrix.</param>
        /// <returns>The converted result.</returns>
        public static double4x4 GetXzyEcefFromXzyEnuMatrix(GeodeticCoordinates origin)
        {
            TrigonometricRatios ratios = new TrigonometricRatios(origin);

            double3 xzyecefPosition = GeodeticToXzyEcef(in ratios, origin.Elevation);

            return new double4x4(
                -ratios.SinLon, ratios.CosLon * ratios.CosLat, -ratios.CosLon * ratios.SinLat, xzyecefPosition.x,
                           0.0,                 ratios.SinLat,                  ratios.CosLat, xzyecefPosition.y,
                 ratios.CosLon, ratios.SinLon * ratios.CosLat, -ratios.SinLon * ratios.SinLat, xzyecefPosition.z,
                           0.0,                           0.0,                            0.0,               1.0
            );
        }

        /// <summary>
        /// Convert a euclidean transformation to a geodetic (latitude/longitude) format.
        /// </summary>
        /// <param name="xzyEcef">The position to convert.</param>
        /// <returns>The converted result.</returns>
        public static GeodeticCoordinates GetGeodeticCoordinates(double3 xzyEcef)
        {
            //see https://hal.science/hal-01704943/file/AccurateEcefConversion-31oct2019.pdf
            double k_CUBICROOT2 = math.pow(2, 1 / 3.0);
            double lat, lon, alt;
            double x = xzyEcef.x;
            double y = xzyEcef.z;
            double z = xzyEcef.y;

            double w2 = x * x + y * y;
            double m = w2 * k_INVA2;
            double n = (z * z) * k_P1ME2 * k_INVA2;
            double p = (m+n - 4 * k_HALFE2 * k_HALFE2) / 6;
            double G = m * n * k_HALFE2 * k_HALFE2;
            double H = 2 * p * p * p + G;
            double C = math.pow(H + G + 2 * math.sqrt(H * G), 1.0 / 3.0) / k_CUBICROOT2;
            double i = -(2 * k_HALFE2 * k_HALFE2 + m + n) / 2;
            double P = p * p;
            double beta = i / 3 - C - P / C;
            double k = k_HALFE2 * k_HALFE2 * (k_HALFE2 * k_HALFE2 - m - n);

            // Compute left part of t
            double t1 = beta * beta - k;
            double t2 = math.sqrt(t1);
            double t3 = t2 - 0.5 * (beta + i);
            double t4 = math.sqrt(t3);
            // Compute right part of t
            double t5 = 0.5 * (beta - i);
            // t5 may accidentally drop just below zero due to numeric turbulence
            // This only occurs at latitudes close to +- 45.3 degrees
            t5 = math.abs(t5);
            double t6 = math.sqrt(t5);
            double t7 = (m < n) ? t6 : -t6;
            // Add left and right parts
            double t = t4 + t7;

            double F = math.pow(t, 4) + 2 * i * t * t + 2 * k_HALFE2 * (m-n) * t + k;
            double dFdt = 4 * math.pow(t, 3) + 4 * i * t + 2 * k_HALFE2 * (m-n);
            double dt = - F / dFdt;
            double u = t + dt + k_HALFE2;
            double v = t + dt - k_HALFE2;
            double w = math.sqrt(w2);

            double dw = w * (1 - 1 / u);
            double dz = z * (1 - k_P1ME2 / v);

            double da = math.sqrt(dw * dw + dz * dz);

            lat = math.atan2(z*u, w*v);
            lon = math.atan2(y, x);
            alt = (u < 1) ? -da : da;

            lat = math.degrees(lat);
            lon = math.degrees(lon);

            return new GeodeticCoordinates(lat, lon, alt);

        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Square(double d)
        {
            return d * d;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Cubic(double d)
        {
            return d * d * d;
        }

        /// <summary>
        /// Convert a euclidean transformation with its rotation to a geodetic (latitude/longitude) format.
        /// </summary>
        /// <param name="position">The translation of the vector to convert.</param>
        /// <param name="rotation">The orientation of the vector to convert.</param>
        /// <returns>The converted result.</returns>
        public static GeodeticTR XzyEcefToGeodetic(double3 position, quaternion rotation)
        {
            GeodeticTR result;

            result.Position = GetGeodeticCoordinates(position);

            double4x4 ecefFromEnu = GetXzyEcefFromXzyEnuMatrix(result.Position);

            quaternion ecefFromEnuRotation = UGMath.ValidTRS(ecefFromEnu)
                ? ecefFromEnu.GetRotation()
                : quaternion.identity;

            quaternion enuFromEcefRotation = math.inverse(ecefFromEnuRotation);
            quaternion enuRotation = math.mul(enuFromEcefRotation, rotation);

            result.EulerAngles = FlipPrincipalAxes(((Quaternion)enuRotation).eulerAngles);

            return result;
        }

        /// <summary>
        /// Validate a coordinate tile is in the valid WGS84 limits.
        /// </summary>
        /// <param name="minLat">Minimum latitude</param>
        /// <param name="maxLat">Maximum latitude</param>
        /// <param name="minLon">Minimum longitude</param>
        /// <param name="maxLon">Maximum longitude</param>
        /// <param name="minEle">Minimum elevation</param>
        /// <param name="maxEle">Maximum elevation</param>
        /// <exception cref="InvalidOperationException">If the coordinates are outside the limits.</exception>
        private static void ValidateCoordinates(double minLat, double maxLat, double minLon, double maxLon, double minEle, double maxEle)
        {
            if (minLon < -180.0 || minLon >= 180.0)
                throw new InvalidOperationException("Region Bounding Volume cannot have west value outside of -PI to PI range");

            if (maxLon <= -180.0 || maxLon > 180.0)
                throw new InvalidOperationException("Region Bounding Volume cannot have east value outside of -PI to PI range");

            if (maxLat < minLat)
                throw new InvalidOperationException("Region Bounding Volume cannot have north value smaller than south value");

            if (maxEle < minEle)
                throw new InvalidOperationException("Region Bounding Volume cannot have max height value smaller than min height");
        }

        /// <summary>
        /// Convert a region as defined by lat, lon and elevation into an axis aligned bounding volume
        /// in ECEF space.
        /// </summary>
        /// <param name="minLat">Minimum latitude, in degrees</param>
        /// <param name="maxLat">Maximum latitude, in degrees</param>
        /// <param name="minLon">Minimum longitude, in degrees</param>
        /// <param name="maxLon">Maximum longitude, in degrees</param>
        /// <param name="minEle">Minimum elevation, in meters</param>
        /// <param name="maxEle">Maximum elevation, in meters</param>
        /// <returns></returns>
        public static unsafe DoubleBounds ConvertRegionBoundingVolume(double minLat, double maxLat, double minLon, double maxLon, double minEle, double maxEle)
        {
            ValidateCoordinates(minLat, maxLat, minLon, maxLon, minEle, maxEle);

            while (maxLon < minLon)
                maxLon += 360.0;

            //
            //  What I'm calling prime latitude is the latitude that is closest
            //      to the equator in the region indicated within this bounding
            //      volume.
            //
            double primeLatitude;
            if (maxLat > 0 && minLat > 0)
                primeLatitude = minLat;
            else if (minLat < 0 && maxLat < 0)
                primeLatitude = maxLat;
            else
                primeLatitude = 0;

            double3* extremes = stackalloc double3[16];
            int extremesLength = 0;

            if (IsMonotonous(minLon, -180.0, maxLon) || IsMonotonous(minLon, 180.0, maxLon))
                extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, -180.0, maxEle);

            if (IsMonotonous(minLon, -90.0, maxLon) || IsMonotonous(minLon, 270.0, maxLon))
                extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, -90.0, maxEle);

            if (IsMonotonous(minLon, 0.0, maxLon) || IsMonotonous(minLon, 360.0, maxLon))
                extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, 0.0, maxEle);

            if (IsMonotonous(minLon, 90.0, maxLon) || IsMonotonous(minLon, 450.0, maxLon))
                extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, 90.0, maxEle);

            extremes[extremesLength++] = GeodeticToXzyEcef(minLat, minLon, minEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(minLat, maxLon, minEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(maxLat, minLon, minEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(maxLat, maxLon, minEle);

            extremes[extremesLength++] = GeodeticToXzyEcef(minLat, minLon, maxEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(minLat, maxLon, maxEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(maxLat, minLon, maxEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(maxLat, maxLon, maxEle);

            extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, minLon, maxEle);
            extremes[extremesLength++] = GeodeticToXzyEcef(primeLatitude, maxLon, maxEle);

            double3 min = new double3(double.MaxValue, double.MaxValue, double.MaxValue);
            double3 max = new double3(double.MinValue, double.MinValue, double.MinValue);

            for (int i = 0; i < extremesLength; i++)
            {
                min = math.min(min, extremes[i]);
                max = math.max(max, extremes[i]);
            }

            return new DoubleBounds(0.5 * (min + max), max - min);
        }

        private static bool IsMonotonous(double a, double b, double c)
        {
            return (a <= b && b <= c);
        }

        private static float3 FlipPrincipalAxes(float3 input)
        {
            return new float3(input.x * -1, input.y, input.z * -1);
        }

        private static void FlipPrincipalAxes(ref quaternion input)
        {
            input.value.x *= -1;
            input.value.z *= -1;
        }
    }
}
