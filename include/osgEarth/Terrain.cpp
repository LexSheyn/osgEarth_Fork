// TEST:
#include "../../include/osgEarth/Terrain.hpp"
//

#include "osgEarth/Terrain"

#include <osgViewer/Viewer>

#define LC "[Terrain] "

namespace osgEarth
{
    ////////////////////////////////////////////////////////////
    ////////// TerrainOperation
    ////////////////////////////////////////////////////////////
    TerrainOperation::TerrainOperation(
        const std::string& name,
        const TileKey&     key,
        osg::Node*         p_tile,
        Terrain*           p_terrain
    )
    : osg::Operation(name, true)
    , m_key     (key)
    , m_tile    (p_tile)
    , m_terrain (p_terrain)
    , m_count   (0)
    , m_delay   (0)
    {}

    const TileKey& TerrainOperation::getKey() const
    {
        return m_key;
    }

    const osg::observer_ptr<osg::Node>& TerrainOperation::getTile() const
    {
        return m_tile;
    }

    const osg::observer_ptr<Terrain>& TerrainOperation::getTerrain() const
    {
        return m_terrain;
    }

    void TerrainOperation::operator()(osg::Object* p_object)
    {
        if (!this->getKeep())
        {
            return;
        }

        if (m_delay > 0)
        {
            m_delay--;

            return;
        }

        m_count++;

        this->execute(p_object);

        this->setKeep(false);
    }

    uint32_t TerrainOperation::getCount() const
    {
        return m_count;
    }

    int32_t TerrainOperation::getDelay() const
    {
        return m_delay;
    }

    void TerrainOperation::setDelay(int32_t delay)
    {
        m_delay = delay;
    }

    ////////////////////////////////////////////////////////////
    ////////// OnTileUpdateOperation
    ////////////////////////////////////////////////////////////
    OnTileUpdateOperation::OnTileUpdateOperation(
        const TileKey& key,
        osg::Node*     p_tile,
        Terrain*       p_terrain
    )
    : TerrainOperation("OnTileUpdateOperation", key, p_tile, p_terrain)
    {}

    void OnTileUpdateOperation::execute(osg::Object* p_object)
    {
        TileKey key = this->getKey();

        osg::ref_ptr<osg::Node> tile;
        osg::ref_ptr<Terrain>   terrain;

        const bool b_keyIsValid    = key.valid();
        const bool b_tileExists    = this->getTile().lock(tile);
        const bool b_terrainExists = this->getTerrain().lock(terrain);

        if (b_keyIsValid && b_tileExists)
        {
            if (b_keyIsValid)
            {
                terrain->onTileUpdated(key, tile.get());
            }
            else
            {
                terrain->onTileUpdated(key, nullptr);
            }
        }
        else
        {
            // Tile expired - let it go.
            OSG_DEBUG_FP << "Tile has expired before notification: " << key.str() << std::endl;
        }
    }

    ////////////////////////////////////////////////////////////
    ////////// Terrain
    ////////////////////////////////////////////////////////////
    bool Terrain::getHeight(
        osg::Node*              p_patch,
        const SpatialReference* p_srs,
        double                  x,
        double                  y,
        double*                 out_heightAboveMSL,
        double*                 out_heightAboveEllipsoid
    ) const
    {
        if (!m_graph.valid() && !p_patch)
        {
            return false;
        }

        // Convert to map coordinates.

        const SpatialReference* p_terrainSRS = this->getSRS();

        if (p_srs && !p_srs->isHorizEquivalentTo(p_terrainSRS))
        {
            p_srs->transform2D(x, y, p_terrainSRS, x, y);
        }

        // Trivially reject a point that lies outside the terrain.

        if (!this->getProfile()->getExtent().contains(x, y))
        {
            return false;
        }

        if (p_srs && p_srs->isGeographic())
        {
            // Pertrub polar latitudes slightly to prevent intersection anomaly at the poles.

            if (p_terrainSRS->isGeographic())
            {
                if (osg::equivalent(y, 90.0))
                {
                    y -= 1e-7;
                }
                else if (osg::equivalent(y, -90.0))
                {
                    y += 1e-7;
                }
            }
        }

        const osg::EllipsoidModel* p_ellipsoidModel = p_terrainSRS->getEllipsoid();

        double radius = osg::minimum(p_ellipsoidModel->getRadiusEquator(), p_ellipsoidModel->getRadiusPolar());

        // Calculate the endpoints for an intersection test.

        osg::Vec3d begin(x, y, radius);
        osg::Vec3d end  (x, y, -radius);

        if (p_terrainSRS->isGeographic())
        {
            const SpatialReference* p_ECEF = p_terrainSRS->getGeocentricSRS();

            p_terrainSRS->transform(begin, p_ECEF, begin);
            p_terrainSRS->transform(end,   p_ECEF, end);
        }

        osgUtil::LineSegmentIntersector* lineSegmentIntersector = new osgUtil::LineSegmentIntersector(begin, end);

        lineSegmentIntersector->setIntersectionLimit(osgUtil::Intersector::LIMIT_NEAREST);

        osgUtil::IntersectionVisitor intersectionVisitor(lineSegmentIntersector);

        if (p_patch)
        {
            p_patch->accept(intersectionVisitor);
        }
        else
        {
            m_graph->accept(intersectionVisitor);
        }

        osgUtil::LineSegmentIntersector::Intersections& results = lineSegmentIntersector->getIntersections();

        if (!results.empty())
        {
            const osgUtil::LineSegmentIntersector::Intersection& firstHit = *results.begin();

            osg::Vec3d hitPoint = firstHit.getWorldIntersectPoint();

            *out_heightAboveMSL = hitPoint.z();

            p_terrainSRS->transformFromWorld(hitPoint, hitPoint, out_heightAboveEllipsoid);

            return true;
        }

        return false;
    }

    bool Terrain::getHeight(
        const SpatialReference* p_srs,
        double                  x,
        double                  y,
        double*                 out_heightAboveMSL,
        double*                 out_heightAboveEllipsoid
    ) const
    {
        return this->getHeight(nullptr, p_srs, x, y, out_heightAboveMSL, out_heightAboveEllipsoid);
    }

    bool Terrain::getWorldCoordsUnderMouse(
        osg::View*  p_view,
        float32_t   mx,
        float32_t   my,
        osg::Vec3d& out_worldCoords
    ) const
    {
        osgViewer::View* p_viewerView = dynamic_cast<osgViewer::View*>(p_view);

        if (!p_viewerView || !m_graph.valid())
        {
            return false;
        }

        float32_t y = 0.0f;
        float32_t x = 0.0f;

        const osg::Camera* p_camera = p_viewerView->getCameraContainingPosition(mx, my, x, y);

        if (!p_camera)
        {
            p_camera = p_viewerView->getCamera();
        }

        // Build a matrix that transforms from the terrain/ world space
        // to either clip of window space, depending on whether or not we
        // have a viewport. Is it even possible to not have a viewport? -gw

        osg::Matrixd matrix;

        // Compensate for any transforms applied between terrain and camera.

        osg::Matrix terrainReferenceFrame = osg::computeLocalToWorld(m_graph->getParentalNodePaths()[0]);

        matrix.postMult(terrainReferenceFrame);
        matrix.postMult(p_camera->getViewMatrix());
        matrix.postMult(p_camera->getProjectionMatrix());

        float64_t zNear = -1.0;
        float64_t zFar  = 1.0;

        if (const osg::Viewport* p_viewport = p_camera->getViewport())
        {
            matrix.postMult(p_viewport->computeWindowMatrix());

            zNear = 0.0;
            zFar  = 1.0;
        }

        osg::Matrixd inverse = osg::Matrixd::inverse(matrix);

        osg::Vec3d beginVertex = osg::Vec3d(x, y, zNear) * inverse;
        osg::Vec3d endVertex   = osg::Vec3d(x, y, zFar)  * inverse;

        osg::ref_ptr<osgUtil::LineSegmentIntersector> lineSegmentIntersector = new osgUtil::LineSegmentIntersector(osgUtil::Intersector::MODEL, beginVertex, endVertex);

        // Limit it to one intersection. We only care about the nearest.

        lineSegmentIntersector->setIntersectionLimit(osgUtil::Intersector::LIMIT_NEAREST);

        osgUtil::IntersectionVisitor intersectionVisitor(lineSegmentIntersector.get());

        m_graph->accept(intersectionVisitor);

        if (lineSegmentIntersector->containsIntersections())
        {
            out_worldCoords = lineSegmentIntersector->getIntersections().begin()->getWorldIntersectPoint();

            return true;
        }

        return false;
    }

    void Terrain::addTerrainCallback(TerrainCallback* p_callback)
    {
        if (p_callback)
        {
            this->removeTerrainCallback(p_callback);

            Threading::ScopedWriteLock exclusiveLock(m_callbackMutex);

            m_callbacks.push_back(p_callback);

            m_callbackCount++;
        }
    }

    void Terrain::removeTerrainCallback(TerrainCallback* p_callback)
    {
        Threading::ScopedWriteLock exclusiveLock(m_callbackMutex);

        for (CallbackList::iterator it = m_callbacks.begin(); it != m_callbacks.end();)
        {
            if (it->get() == p_callback)
            {
                it = m_callbacks.erase(it);

                m_callbackCount--;
            }
            else
            {
                it++;
            }
        }
    }

    void Terrain::accept(osg::NodeVisitor& nodeVisitor)
    {
        m_graph->accept(nodeVisitor);
    }

    void Terrain::notifyTileUpdated(
        const TileKey& key,
        osg::Node*     p_tile
    )
    {
        if (!p_tile)
        {
            OSG_WARN << LC << "node is null!" << std::endl;
        }

        if (m_callbackCount > 0)
        {
            if (!key.valid())
            {
                OSG_WARN << LC << "key is invalid!" << std::endl;
            }

            OnTileUpdateOperation* p_onTileUpdateOperation = new OnTileUpdateOperation(key, p_tile, this);

            m_updateQueue->add(p_onTileUpdateOperation);
        }
    }

    void Terrain::notifyTilesRemoved(const std::vector<TileKey>& keys)
    {
        // Not implemented.
    }

    void Terrain::notifyMapElevationChanged()
    {
        if (m_callbackCount > 0)
        {
            OnTileUpdateOperation* p_onTileUpdateOperation = new OnTileUpdateOperation(TileKey::INVALID, nullptr, this);

            // Let the terrain update before applying this.

            p_onTileUpdateOperation->setDelay(1);

            m_updateQueue->add(p_onTileUpdateOperation);
        }
    }

    Terrain::Terrain(
        osg::Node*     p_graph,
        const Profile* p_mapProfile
    )
    : m_callbackMutex (OE_MUTEX_NAME)
    , m_callbackCount (0)
    , m_graph         (p_graph)
    , m_profile       (p_mapProfile)
    , m_updateQueue   (new osg::OperationQueue())
    {}

    void Terrain::update()
    {
        m_updateQueue->runOperations();
    }

    void Terrain::onTileUpdated(
        const TileKey& key,
        osg::Node*     p_tile
    )
    {
        Threading::ScopedReadLock sharedLock(m_callbackMutex);

        for (CallbackList::iterator it = m_callbacks.begin(); it != m_callbacks.end();)
        {
            TerrainCallbackContext terrainCallbackContext(this);

            it->get()->onTileUpdate(key, p_tile, terrainCallbackContext);

            // If the callback has "remove" flag set - discard the callback.

            if (terrainCallbackContext.markedForRemoval())
            {
                it = m_callbacks.erase(it);
            }
            else
            {
                it++;
            }
        }
    }

    void Terrain::onTilesRemoved(const std::vector<TileKey>& keys)
    {
        // TO DO
    }

    void Terrain::onMapElevationChanged()
    {
        // TO DO
    }
} // namespace osgEarth