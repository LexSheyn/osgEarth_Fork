#ifndef OSGEARTH_SPATIAL_REFERENCE_H
#define OSGEARTH_SPATIAL_REFERENCE_H 1

#include <osgEarth/Common>
#include <osgEarth/Units>
#include <osgEarth/VerticalDatum>
#include <osgEarth/Threading>
#include <osgEarth/Containers>

#include <osg/CoordinateSystemNode>
#include <osg/Vec3>

#include <unordered_map>
#include <cstdint>

typedef float  float32_t;
typedef double float64_t;

namespace osgEarth
{
    // Definitions for the mercator extent.

    const float64_t MERC_MINX   = -20037508.34278925;
    const float64_t MERC_MINY   = -20037508.34278925;
    const float64_t MERC_MAXX   =  20037508.34278925;
    const float64_t MERC_MAXY   =  20037508.34278925;
    const float64_t MERC_WIDTH  =  MERC_MAXX - MERC_MINX;
    const float64_t MERC_HEIGHT =  MERC_MAXY - MERC_MINY;

    // SpatialReference holds information describing the reference ellipsoid/datum
    // and the projection of geospatial data.
    class OSGEARTH_EXPORT SpatialReference : public osg::Referenced
    {
    public:
        // Creates an SRS from two initialization strings. The first one for the horizontal datum and
        // the second one for the vertical datum. If you omit the vertical datum, it will default to
        // the geodetic datum for the ellipsoid.
        static SpatialReference* create(const std::string& horizontalDatum, const std::string& verticalDatum);

        // Creates and SRS by cloning a pre-existing OGR spatial reference handle.
        // The new SRS owns the cloned handle, and the caller retains responsibility
        // for managing the original handle.
        static SpatialReference createFromHandle(void* ogrHandle);

    public: // Point transformations.

    // Transforms a single point from this SRS to another SRS.
    // @returns
    //      true if the transformation was successfull, otherwise false.
    virtual bool transformPoint(
        const osg::Vec3d&       point,
        const SpatialReference* to_SRS,
        osg::Vec3d*             out_point
    ) const;

    // Transforms a collection of poits from this SRS to another SRS.
    // @returns
    //      trus if all transformations were successfull, otherwise false.
    virtual bool transformPoints(
        std::vector<osg::Vec3d>& points,
        const SpatialReference*  to_SRS
    ) const;

    // Transforms a 2D point directly. (Convenience function)
    bool transformPoint2D(
        float64_t               x,
        float64_t               y,
        const SpatialReference* to_SRS,
        float64_t*              out_x,
        float64_t*              out_y
    ) const;

    public: // Units conversions.

        // Transforms a distance from the base units of this SRS to the base units of
        // another one. If one of the SRS's is geographic (i.e. has angular units), the
        // conversion will assume that the corresponding distance is measured at the
        // equator.
        float64_t convertUnits(
            float64_t               distance,
            const SpatialReference* to_SRS,
            float64_t               latitude
        ) const;

        static float64_t transformUnits(
            const Distance&         distance,
            const SpatialReference* to_SRS,
            float64_t               latitude
        );

    public: // World transformations.

        // Transforms a point from this SRS into "world" coordinates. This normalizes
        // the Z coordinate (according to the vertical datum) and converts to geocentric
        // coordinates if necessary.
        bool transformToWorldCoords(
            const osg::Vec3d& point,
            osg::Vec3*        out_worldPoint
        ) const;

        // Transforms a point from the "world" coordinate system into this spatial
        // reference.
        // @param worldPoint
        //      Point in "world" coordinate system.
        // @param out_localPoint
        //      Output point in local (SRS) coordinates.
        // @param worldIsGeocentric
        //      Whether or not the incoming world coordinates are geocentric.
        // @param out_geodeticZ
        //      (Optional) Outputs the geodetic (HAE) Z if applicable.
        bool transformFromWorldCoords(
            const osg::Vec3d& worldPoint,
            osg::Vec3d*       out_localPoint,
            float64_t*        out_geodeticZ
        ) const;

    public: // Extent transformations.

        // Transforms a spatial extent to another SRS. The transformed extent will
        // actually be minimum bounding axis-aligned rectangle that would hold
        // the source extent.
        virtual bool transformExtentToMBR(
            const SpatialReference* to_SRS,
            float64_t*              out_x_min,
            float64_t*              out_y_min,
            float64_t*              out_x_max,
            float64_t*              out_y_max
        ) const;

        virtual bool transformExtentPoints(
            const SpatialReference* to_SRS,
            float64_t               x_min,
            float64_t               y_min,
            float64_t               x_max,
            float64_t               y_max,
            float64_t*              out_x,
            float64_t*              out_y,
            uint32_t                countX,
            uint32_t                countY
        ) const;
    };

    //public: // Properties.

    // TO DO:
    // File: SpatialReference
    // Line: 165
}

#endif // OSGEARTH_SPATIAL_REFERENCE_H