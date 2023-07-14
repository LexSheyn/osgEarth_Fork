#include "osgEarth/Terrain"

#include <osgViewer/Viewer>

namespace osgEarth
{
    Terrain::OnTileUpdateOperation::OnTileUpdateOperation(
        const TileKey& key,
        osg::Node*     p_node,
        Terrain*       p_terrain,
    )
    : osg::Operation()
    , m_key     (key)
    , m_node    (p_node)
    , m_terrain (p_terrain)
    , m_count   (0)
    , m_delay   (0)
    {}

    void Terrain::OnTileUpdateOperation::operator()(osg::Object* p_object)
    {
        if (!getKeep())
        {
            return;
        }
        
        if (m_delay > 0)
        {
            m_delay--;

            return;
        }

        m_count++;

        osg::ref_ptr<Terrain>   terrain;
        osg::ref_ptr<osg::Node> node;

        if (m_terrain.lock(terrain) && (!m_key.valid()))
        {
            m_terrain->notifyTileUpdate(m_key, nullptr);
        }
    }
}