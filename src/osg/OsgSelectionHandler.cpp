#include "OsgSelectionHandler.hpp"
#include <osgUtil/LineSegmentIntersector>
#include <osg/MatrixTransform>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/PolygonMode>
#include <osg/ValueObject>
#include <iostream>

void OsgSelectionHandler::HandleSelection(osg::Camera* cam, int x, int y) {

    osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector =
            new osgUtil::LineSegmentIntersector(
                osgUtil::Intersector::WINDOW,
                static_cast<double>(x),
                static_cast<double>(y));

    osgUtil::IntersectionVisitor iv(intersector.get());
    iv.setTraversalMask(~0x1);
    cam->accept(iv);

    if(intersector->containsIntersections())
    {
        const osgUtil::LineSegmentIntersector::Intersection& result =
                *(intersector->getIntersections().begin());

        osg::BoundingBox bb = result.drawable->getBoundingBox();
        auto selection_box_id = process_selections(result.drawable->getParent(0));

        if(selection_box_id != -1)
        {
            osg::MatrixTransform* active_selection_box = dynamic_cast<osg::MatrixTransform*>(m_selection_boxes->getChild(selection_box_id));
            if(active_selection_box)
            {
                osg::Vec3 worldCenter = bb.center() * osg::computeLocalToWorld(result.nodePath);
                active_selection_box->setMatrix(osg::Matrix::scale(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin()) * osg::Matrix::translate(worldCenter));
            }
            else
            {
                std::cout << "ERROR: Active selection box error!" << std::endl;
            }
        }
    }
}

osg::Switch* OsgSelectionHandler::GetOrCreateSelectionBoxSwitch() {

    if(!m_selection_boxes)
        m_selection_boxes = new osg::Switch;
    return m_selection_boxes.get();
}


int OsgSelectionHandler::process_selections(osg::Node* node) {

    bool selection = false;
    if(node->getUserValue("Selection", selection)) {
        // invert the selection
        selection = !selection;
        int selection_box_id = -1;
        if(selection) {
            // select the node
            node->setUserValue("Selection", selection);
            selection_box_id = get_or_create_selection_box();
            node->setUserValue("Selection_Box_Id", selection_box_id);
            m_selection_boxes->setValue(static_cast<unsigned int>(selection_box_id), selection);
            m_selection_boxes->getChild(static_cast<unsigned int>(selection_box_id))->setUserValue("Free", !selection);
            return selection_box_id;
        }
        else {
            // deselect the node
            if(node->getUserValue("Selection_Box_Id", selection_box_id)) {
                node->setUserValue("Selection", selection);
                node->setUserValue("Selection_Box_Id", -1);
                m_selection_boxes->setValue(static_cast<unsigned int>(selection_box_id), selection);
                m_selection_boxes->getChild(static_cast<unsigned int>(selection_box_id))->setUserValue("Free", !selection);
            }
            else {
                std::cout << "ERROR: This node does not have Selection_Box_Id" << std::endl;
            }
        }
    }
    else {
        if(node->getNumParents() != 0) {
            return process_selections(node->getParent(0));
        }
    }
    return -1;
}

unsigned int OsgSelectionHandler::get_or_create_selection_box() {

    for(unsigned int i = 0; i < m_selection_boxes->getNumChildren(); ++i)
        if(is_free(m_selection_boxes->getChild(i)))
            return i;

    return create_selection_box(1);
}

bool OsgSelectionHandler::is_free(osg::Node* node) {

    bool free = false;
    if(node->getUserValue("Free", free))
        return free;
    else
        std::cout << "ERROR: node is not a selection box" << std::endl;
    return false;
}

unsigned int OsgSelectionHandler::create_selection_box(int num) {

    for(int i = 0; i < num; ++i) {
        osg::ref_ptr<osg::Geode> geode = new osg::Geode;
        osg::ShapeDrawable* shape = new osg::ShapeDrawable(new osg::Box(osg::Vec3(), 2.0f));
        shape->setColor(osg::Vec4d(1.0f, 0.0f, 0.0f, 1.0f));
        geode->addDrawable(shape);
        osg::ref_ptr<osg::MatrixTransform> selection_box = new osg::MatrixTransform;
        selection_box->setNodeMask(0x1);
        selection_box->addChild(geode.get());
        osg::StateSet* ss = selection_box->getOrCreateStateSet();
        ss->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
        ss->setAttributeAndModes(new osg::PolygonMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE));
        selection_box->setUserValue("Free", true);
        m_selection_boxes->addChild(selection_box.release(), false);
    }

    return m_selection_boxes->getNumChildren() - 1;
}
