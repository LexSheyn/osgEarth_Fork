#ifndef OSGEARTH_TERRAIN_H
#define OSGEARTH_TERRAIN_H 1

#include "osgEarth/Common"
#include "osgEarth/TileKey"
#include "osgEarth/Threading"

#include <osg/OperationThread>
#include <osg/View>

#include <iostream>
#include <cstdint>
#include <list>
#include <atomic>

typedef float  float32_t;
typedef double float64_t;

namespace osgEarth
{
    class Terrain;
    class SpatialReference;

    // This ibject is passed to terrain callbacks to provide context information
    // and a feedback interface to terrain callback consumers.
    class TerrainCallbackContext
    {
    public:
        const Terrain* getTerrain() const { return mp_terrain; }

        // Indicate that the callback should be removed immediately.
        void remove() { mb_remove = true; }

        // Whether the used called remove().
        bool markedForRemoval() const { return mb_remove; }

    public:
        TerrainCallbackContext(Terrain* p_terrain)
        : mb_remove  (false)
        , mp_terrain (p_terrain)
        {}

        virtual ~TerrainCallbackContext() {}

    protected:
        bool     mb_remove;
        Terrain* mp_terrain;

        friend class Terrain;
    };

    // Callback that you can register with the Terrain in order to receive
    // update messages about the terrain scene graph.
    class TerrainCallback : public osg::Referenced
    {
    public:
        // A tile was added to the terrain graph.
        // @param key
        //      Tile key of the new tile, including the geographic extents.
        // @param graph
        //      Scene graph that can be used to intersect the tile.
        //      It may include more than just new tile.
        // @param context
        //      Contextual information about the callback.
        virtual void onTileUpdate(
            const TileKey&          key,
            osg::Node*              p_graph,
            TerrainCallbackContext& context
        )
        {
            // Empty by default.
            // Was calling empty onTileAdded(key, graph, context) for backward compatibility.
        };

        virtual ~TerrainCallback() {}
    };

    // Convenience adapter for hooking a callback into a class. The
    // class must be derived from osg::Referenced, When the calss
    // destructs, the callback will automatically remove itself from
    // the Terrain object.
    template<typename T>
    class TerrainCallbackAdapter : public TerrainCallback
    {
    public:
        TerrainCallbackAdapter(T* t)
        : m_t(t)
        {}

        virtual void onTileUpdate(
            const TileKey&          key,
            osg::Node*              p_tile,
            TerrainCallbackContext& context
        ) override
        {
            osg::ref_ptr<T> t;

            if (m_t.lock(t))
            {
                t->onTileUpdate(key, p_tile, context);
            }
            else
            {
                context.remove();
            }
        }

    protected:
        virtual ~TerrainCallbackAdapter() {}

        osg::observer_ptr<T> m_t;
    };

    // Interface for an object that can resolve the terrain elevation
    // at a given map location.
    class ITerrainResolver
    {
    public:
        // Intersects the terrain at the location x, y and returns the height data.
        // @param p_patch
        //      Subgraph patch.
        // @param srs
        //      Spatial reference system of (x, y) coordinates.
        // @param x
        //      Coordinate x at which to query the height.
        // @param y
        //      Coordinate y at which to query the height.
        // @param out_heightAboveMSL
        //      The resulting relative height goes here. The height is relative to MSL
        //      (mean sea level) as expressed by the map's vertical datum.
        // @param out_heightAboveEllipsoid
        //      The resulting geodetic height goer here. The height is relative to the
        //      geodetic ellipsoid expressed by the map's SRS.
        virtual bool getHeight(
            osg::Node*              p_patch,
            const SpatialReference* p_srs,
            double                  x,
            double                  y,
            double*                 out_heightAboveMSL,
            double*                 out_heightAboveEllipsoid
        ) const = 0;

        // Same as above, but with no subgraph patch provided.
        virtual bool getHeight(
            const SpatialReference* p_srs,
            double                  x,
            double                  y,
            double*                 out_heightAboveMSL,
            double*                 out_heightAboveEllipsoid
        ) const = 0;

        virtual ~ITerrainResolver() {}
    };

    // TerrainHeightProvider is replaced by TerrainResolver (now ITerrailResolver).

    // Base class for the Terrain related operations.
    class TerrainOperation : public osg::Operation
    {
    protected:
        TerrainOperation(
            const std::string& name,
            const TileKey&     key,
            osg::Node*         p_tile,
            Terrain*           terrain
        );

        // Function gets called from within the operator() when delay has expired.
        virtual void execute(osg::Object* p_object) = 0;

        // @returns
        //      Tile key.
        const TileKey& getKey() const;

        // @returns
        //      Tile node observer pointer.
        const osg::observer_ptr<osg::Node>& getTile() const;

        // @returns
        //      Terrain observer pointer.
        const osg::observer_ptr<Terrain>& getTerrain() const;

    public:
        // Calls execute() function only when delay has expired.
        // Otherwise decrements the delay.
        virtual void operator()(osg::Object* p_object) override final;

        // @returns
        //      Number of times this operation was called already.
        uint32_t getCount() const;

        // @returns
        //      Remaning number of update calls before this operation will be executed.
        int32_t getDelay() const;

        // @param
        //      Number of update calls before this operation will be executed.
        void setDelay(int32_t delay);

    private:
        TileKey                      m_key;
        osg::observer_ptr<osg::Node> m_tile;
        osg::observer_ptr<Terrain>   m_terrain;
        uint32_t                     m_count;
        int32_t                      m_delay;
    };

    // Performs sequential tile update callbacks execution.
    class OnTileUpdateOperation : public TerrainOperation
    {
    public:
        OnTileUpdateOperation(
            const TileKey& key,
            osg::Node*     p_tile,
            Terrain*       terrain
        );

        virtual void execute(osg::Object* p_object) override;
    };

    // Services for interacting with the live terrain graph. This differs from
    // the Map model. Map represents the parametric data backing the terrain,
    // while Terrain represents the actual geometry in memory.
    //
    // All returned map coordinate values are in the units conveyed in the
    // spatial reference at getSRS().
    class OSGEARTH_EXPORT Terrain : public osg::Referenced, public ITerrainResolver
    {
    public:
        // Gets the profile of the map with which this terrain is associated.
        const Profile* getProfile() const { return m_profile.get(); }

        // Gets the spatial reference of the map with which this terrain is assosiated.
        const SpatialReference* getSRS() const { return m_profile->getSRS(); }

    public: // ITerrainResolver interface

        // Intersects the terrain at the location x, y and returns the height data.
        // @param p_patch
        //      Subgraph patch.
        // @param srs
        //      Spatial reference system of (x, y) coordinates.
        // @param x
        //      Coordinate x at which to query the height.
        // @param y
        //      Coordinate y at which to query the height.
        // @param out_heightAboveMSL
        //      The resulting relative height goes here. The height is relative to MSL
        //      (mean sea level) as expressed by the map's vertical datum.
        // @param out_heightAboveEllipsoid
        //      The resulting geodetic height goer here. The height is relative to the
        //      geodetic ellipsoid expressed by the map's SRS.
        virtual bool getHeight(
            osg::Node*              p_patch,
            const SpatialReference* p_srs,
            double                  x,
            double                  y,
            double*                 out_heightAboveMSL,
            double*                 out_heightAboveEllipsoid
        ) const override;

        // Same as above, but with no subgraph patch provided.
        virtual bool getHeight(
            const SpatialReference* p_srs,
            double                  x,
            double                  y,
            double*                 out_heightAboveMSL,
            double*                 out_heightAboveEllipsoid
        ) const override;

    public:

        // Gives the world coordinates under the mouse.
        // @oaram view
        //      View in which to do the query.
        // @param mx
        //      Mouse x coordinate.
        // @param my
        //      Mouse y coordinate.
        // @param out_coords
        //      Stores the world coordinates under the mouse (when returning true).
        // @returns
        //      true if a mouse is over the Terrain, false otherwise.
        bool getWorldCoordsUnderMouse(
            osg::View*  p_view,
            float32_t   mx,
            float32_t   my,
            osg::Vec3d& out_worldCoords
        ) const;

    public:

        // Adds a terrain callback.
        // @param callback
        //      Terrain callback to add. This will get called wheneven tile data changes in
        //      the active terrain graph.
        void addTerrainCallback(TerrainCallback* p_callback);

        // Removes a terrain callbac,
        void removeTerrainCallback(TerrainCallback* p_callback);

    public:
        // Accept a node visitor on the terrains's scene graph.
        void accept(osg::NodeVisitor& nodeVisitor);

        // Access the raw terrain graph.
        osg::Node* getGraph() const { return m_graph.get(); }

        // Queues the onTileUpdate callback (internal).
        void notifyTileUpdated(
            const TileKey& key,
            osg::Node*     p_tile
        );

        // Not implemented.
        // Queues the onTileRemoved callback (internal).
        void notifyTilesRemoved(const std::vector<TileKey>& keys);

        // Currently queues onTileUpdate with 1 cycle delay.
        // Supposed to queue the onMapElevationChanged (internal).
        void notifyMapElevationChanged();

        virtual ~Terrain() {}

    private:
        Terrain(
            osg::Node*     p_graph,
            const Profile* p_mapProfile
        );

        // Update traversal.
        void update();

        // Sequentially calls onTileUpdate for every tile.
        void onTileUpdated(
            const TileKey& key,
            osg::Node*     p_tile
        );

        // Not implemented.
        void onTilesRemoved(const std::vector<TileKey>& keys);

        // Not implemented.
        void onMapElevationChanged();

        friend class TerrainEngineNode;
        friend class OnTileUpdateOperation;
        friend class OnTileRemoveOperation;
        friend class OnMapElevationChangeOperation;

        typedef std::list<osg::ref_ptr<TerrainCallback>> CallbackList;

        CallbackList                      m_callbacks;
        Threading::ReadWriteMutex         m_callbackMutex;
        std::atomic<int32_t>              m_callbackCount; // Separate size tracker for MT size check without a lock.
        osg::observer_ptr<osg::Node>      m_graph;
        osg::ref_ptr<const Profile>       m_profile;
        osg::ref_ptr<osg::OperationQueue> m_updateQueue;
    };
} // namespace osgEarth

#endif // OSGEARTH_TERRAIN_H
