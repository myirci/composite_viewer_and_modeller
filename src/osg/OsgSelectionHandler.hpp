#ifndef OSG_SELECTION_HANDLER_HPP
#define OSG_SELECTION_HANDLER_HPP

#include <osg/Switch>
#include <osg/Camera>

class OsgSelectionHandler {
public:
    OsgSelectionHandler() { }
    osg::Switch* GetOrCreateSelectionBoxSwitch();
    void HandleSelection(osg::Camera* cam, int x, int y);
private:
    osg::ref_ptr<osg::Switch> m_selection_boxes;
    unsigned int create_selection_box(int num);
    int process_selections(osg::Node* node);
    inline bool is_free(osg::Node* node);
    unsigned int get_or_create_selection_box();
};

#endif // OSG_SELECTION_HANDLER_HPP
