'''
utils.py

snapper v 0.0.1
Copyright 2023 Jason Mak

This file is part of snapper.

snapper is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

snapper is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with snapper. If not, see <https://www.gnu.org/licenses/>.

'''

import pya
import math
import snapper.ports as ports

from functools import singledispatch

def get_version():
  '''Return the version as integers
'''
  version = pya.__version__
  major_str, minor_str, micro_str = version.split(".")
  return int(major_str), int(minor_str), int(micro_str)
  
  
def get_layer_info_dict():
  '''
  Returns current layout's <pya.LayerInfo> in <dict> format

  Arguments
  --------

  Returns
  -------
  :layer_info_dict: <dict> of pairing <int> layer index and corresponding <pya.LayerInfo> object

  '''
  lv = pya.Application.instance().main_window().current_view()
  layout = lv.active_cellview().layout()

  layer_info_dict = {}
  for layer in lv.each_layer():
    layer_index = layer.layer_index()
    layer_info_dict[layer_index] = layout.get_info(layer_index)
  return layer_info_dict

def grid(dbu = True):
  '''
  Returns current grid value

  Arguments
  --------

  kwargs
  :dbu: <bool> Flag to return grid size in dbu of the layout

  Returns
  -------

  if :dbu:
  <int> Grid size in dbu

  else:
  <float> Grid size in measurement units (e.g. microns)


  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()

  g_str = pya.Application.instance().get_config('edit-grid')
  g_float = ly.dbu
  
  if any([i.isdigit() for i in g_str]):
    g_float = float(g_str)
    
  if not dbu:
    return g_float
  else:
    return int(round(g_float/ly.dbu))


###
# View, Navigation, Position
###

def zoom(w):
  '''
  Sets zoom to size :w:, retaining the current center of the view

  Arguments
  --------
  :w: <float> Width (in measurement units) of the new view

  GUI Action
  -------
  Sets the zoom box to a square of size 2w centered at the current location

  '''

  lv = pya.Application.instance().main_window().current_view()
  c = lv.box().center()
  new_box = pya.DBox(c.x - w, c.y - w, c.x + w, c.y + w)
  lv.zoom_box(new_box)


def zoom_cursor(w):
  '''
  Sets zoom to size :w:, centered at the current cursor location

  Arguments
  --------
  :w: <float> Width (in measurement units) of the new view

  GUI Action
  -------
  Sets the zoom box to a square of size 2w centered at the current cursor location

  '''
  lv = pya.Application.instance().main_window().current_view()
  p_c = cursor_position(in_dbu=False, relative=False, on_grid=False)
  # c = lv.box().center()
  new_box = pya.DBox(p_c.x - w, p_c.y - w, p_c.x + w, p_c.y + w)
  lv.zoom_box(new_box)


def orientation(u, v, relative=True, as_axis_direction=False):
  '''
  Converts pixel coordinates/axis into view coordinates/axis

  Arguments
  ---------
  :u: <float> horizontal pixel coordinate
  :v: <float> vertical pixel coordinate

  kwargs
  :relative: <bool> Flag for relative coordinate system (if we are descended into a cell instance)
  :as_axis_direction: <bool> Flag to return as a string tag (e.g. "x", "y") and sign (1, -1)


  '''
  ## Direction of viewport
  # u, v = 0,0: bottom left
  # u increasing to right
  # v increasing to top
  # Right: 1, 0
  # Left: -1, 0
  # Up: 0,1
  # Down: 0, -1
  lv = pya.Application.instance().main_window().current_view()
  global_trans = pya.ICplxTrans().from_s(lv.get_config("global-trans"))

  x_g, y_g = rotate_axis(u, v, global_trans)
  if not relative:
    x_r, y_r = x_g, y_g

  cv = lv.active_cellview()
  ctx_trans = pya.ICplxTrans()
  # mag, rot, mirrx, x, y)
  # the rotation and mirroring happens about the origin of the new instance.
  # need to account for that after

  x_p = int(x_g)
  y_p = int(y_g)

  for ie in cv.context_path:  # the current instance path of the active cellview
    ie_trans = ie.specific_cplx_trans()
    ctx_trans = ctx_trans * ie_trans
    # ie_trans.invert()
    # print("trans", ie_trans.to_s()) #prints the string
  x_p, y_p = rotate_axis(x_p, y_p, ctx_trans)
  x_r, y_r = x_p, y_p

  if not as_axis_direction:
    return x_r, y_r

  # else, present in "x/y" and direction
  if abs(x_r) > abs(y_r):
    axis = "x"
    direction = round(x_r / abs(x_r))
  else:
    axis = "y"
    direction = round(y_r / abs(y_r))
  return axis, direction


def filter_qt_children(qt_object, child_type):
  return [a for a in qt_object.children() if type(a) == child_type]

def _viewport_u_v():
  '''Cursor position in the layout viewport
  Implementation from KLayout changed after version 0.26
  '''
  maj, minor, micro = get_version()
  
  if minor<27:
    lv = pya.Application.instance().main_window().current_view()
    cursor = lv.cursor
  
    # current position in global
    pos_global = cursor.pos
  
    # relative to LayoutView
    pos_local = lv.mapFromGlobal(pos_global)
  
    # convenience function to get viewport height and width
    # alternately, use lv.visibleRegion().rects[0].width and lv.visibleRegion().rects[0].height
    pix_height = lv.viewport_height()
    pix_width = lv.viewport_width()
  
    y_pix_offset = 1  # some kind of buffer that was throwing values off

  else:
    #we have to look for the viewport QWidget because layoutview is no longer inherited
    mw = pya.Application.instance().main_window()
    qstackedwidget = filter_qt_children(mw, pya.QStackedWidget_Native)[0] #This is the LayoutView tab object
    
    qframes = filter_qt_children(qstackedwidget, pya.QFrame_Native)
    #the qframe containing the view widget will have the same width and height as the stackedwidget
    #qframe = [a for a in qframes if (a.width == qstackedwidget.width) and (a.height == qstackedwidget.height)][0] #This is the viewport region corresponding to lv.box()
    
    #qframe.children() includes the tab bar, the close button and the view
    view_widget_list = []
    for qf in qframes:
      view_widget_list += filter_qt_children(qf, pya.QWidget_Native)
    view_widget = view_widget_list[0]
      
    
    pos_global = view_widget.cursor.pos
    
    # relative to the view
    pos_local = view_widget.mapFromGlobal(pos_global)
    
    pix_height = view_widget.height
    pix_width = view_widget.width
    y_pix_offset = 1


  # (u,v) = (0,0) at top left
  # v increases top to bottom
  # u increases left to right    
  u = (pos_local.x) / pix_width  # normalized display horizontal coordinate
  v = (pos_local.y + y_pix_offset) / pix_height  # normalized display vertical coordinate
  
  return u, v


def cursor_position(in_dbu=False, relative=True, on_grid=False):
  '''
  Returns current cursor position

  Arguments
  --------
  kwargs
  :in_dbu: <bool> Flag to return in dbu
  :relative: <bool> Flag to return coordinate in of the current cell descended into
  :on_grid: <bool> Flag to return cursor position rounded to the closest grid point

  Returns
  -------
  :cursor_point:
  if :in_dbu:
    <pya.Point>
  else:
    <pya.DPoint>

  if :relative:
    Cursor coordinate in coordinate system of the instance currently descended into
  else
    Global cursor coordinate

  if :on_grid:
    Cursor position closest to the grid point defined by utils.grid()
  else:
    Pixel accurate cursor position
  '''
  lv = pya.Application.instance().main_window().current_view()
  # relative: coordinate system in current cellview instance. Not relative is global coordinate.

  # current global orientation
  global_trans = pya.ICplxTrans().from_s(lv.get_config("global-trans"))

  if on_grid:
    g_dbu = grid()

    # ly = lv.active_cellview().layout()
    # g_str = pya.Application.instance().get_config('edit-grid')
    # g_float = ly.dbu
    # if any([i.isdigit() for i in g_str]):
    #  g_float = float(g_str)
    # g_dbu = int(round(g_float/ly.dbu))

  '''
  cursor = lv.cursor

  # current position in global
  pos_global = cursor.pos

  # relative to LayoutView
  pos_local = lv.mapFromGlobal(pos_global)

  # convenience function to get viewport height and width
  # alternately, use lv.visibleRegion().rects[0].width and lv.visibleRegion().rects[0].height
  pix_height = lv.viewport_height()
  pix_width = lv.viewport_width()

  y_pix_offset = 1  # some kind of buffer that was throwing values off

  # (u,v) = (0,0) at top left
  # v increases top to bottom
  # u increases left to right

  u = (pos_local.x) / pix_width  # normalized display horizontal coordinate
  v = (pos_local.y + 1) / pix_height  # normalized display vertical coordinate
  '''
  
  
  u,v = _viewport_u_v()

  # to convert between global and local pixel positions, refer to Qt 'mapFrom', 'mapFromGlobal', 'mapFromParent', 'mapTo', 'mapToGlobal', 'mapToParent'
  lv_box = lv.box()  # this gives us the properties of the box in microns
  lv_y_span = lv_box.height()
  lv_x_span = lv_box.width()

  # the layoutview box gets rotated relative to what is displayed.
  lv_y_max = lv_box.top
  lv_y_min = lv_box.bottom
  lv_x_max = lv_box.right
  lv_x_min = lv_box.left

  x_micron = 0
  y_micron = 0

  ## case for all 8 transforms
  if global_trans.is_mirror() is False:
    if global_trans.angle == 0:  # r0
      # x = u
      # y = -v

      x_micron = lv_x_min + u * lv_x_span  # top left coordinate + displacement
      y_micron = lv_y_max - v * lv_y_span

      # x_micron = micron_box.left + pos_local.x*micron_width/pix_width
      # y_micron = micron_box.top - (pos_local.y + y_pix_offset)*micron_height/pix_height
      # print(pos_local.x, pos_local.y, x_micron, y_micron)
    elif global_trans.angle == 90:  # r90
      # x = -v
      # y = -u

      x_micron = lv_x_max - v * lv_x_span
      y_micron = lv_y_max - u * lv_y_span

    elif global_trans.angle == 180:  # r180
      # x = -u
      # y = v

      x_micron = lv_x_max - u * lv_x_span
      y_micron = lv_y_min + v * lv_y_span

    elif global_trans.angle == 270:  # r270
      # x = v
      # y = u

      x_micron = lv_x_min + v * lv_x_span
      y_micron = lv_y_min + u * lv_y_span

  else:  # global_trans.is_mirror() is True
    if global_trans.angle == 0:  # m0
      # x = u
      # y = v

      x_micron = lv_x_min + u * lv_x_span  # top left coordinate + displacement
      y_micron = lv_y_min + v * lv_y_span

    elif global_trans.angle == 90:  # m45
      # x = -v
      # y = u

      x_micron = lv_x_max - v * lv_x_span  # top left coordinate + displacement
      y_micron = lv_y_min + u * lv_y_span

    elif global_trans.angle == 180:  # m90
      # x = -u
      # y = -v

      x_micron = lv_x_max - u * lv_x_span  # top left coordinate + displacement
      y_micron = lv_y_max - v * lv_y_span

    elif global_trans.angle == 270:  # m135
      # x = v
      # y = -u

      x_micron = lv_x_min + v * lv_x_span  # top left coordinate + displacement
      y_micron = lv_y_max - u * lv_y_span

  # round to dbu
  ly = lv.active_cellview().layout()
  dbu = ly.dbu
  x_g, y_g = int(round(x_micron / dbu)), int(round(y_micron / dbu))

  if on_grid:  # if grid, round the cursor position to the nearest grid point
    x_g = int(g_dbu * round(x_g / g_dbu))
    y_g = int(g_dbu * round(y_g / g_dbu))

  if not relative:  # if considering only global coordinate system, return values here
    if in_dbu:
      return pya.Point(x_g, y_g)
    else:
      return pya.DPoint(x_g * dbu, y_g * dbu)  # convert back to micron

  # x_g, y_g in dbu at this point

  # else, considering local coordinate system
  cv = lv.active_cellview()
  x_r = int(x_g)  # relative coordinates, in dbu
  y_r = int(y_g)

  # for ie in cv.context_path: #get the path for instances being traversed
  #  x_r, y_r = _utils.coordinate_system_transform(x_r, y_r, ie.specific_cplx_trans()) #get the icplxtrans from the instance, covert the current coordinates into the next coordinate system
  x_r, y_r = relative_coordinates(x_r, y_r, cv.context_path)  # accumulate transforms over the path for instances being traversed

  x_r, y_r = int(round(x_r)), int(round(y_r))

  # return values
  if in_dbu:
    return pya.Point(x_r, y_r)
  else:
    return pya.DPoint(x_r * dbu, y_r * dbu)

###
# Processing selected objects
###

def get_inst_poly_dict(sel_obj):
  '''
    Extracts polygons of a selected object known to be an instance into a <dict>

    Arguments
    --------

    :sel_obj: <pya.ObjectInstPath> Selected object

    Returns
    -------
    :inst_poly_dict: <dict> of keys <int> layer index with values <list> of <pya.Polygon> for the associated layer
    '''
  lv = pya.Application.instance().main_window().current_view()
  layout = lv.active_cellview().layout()

  sel_inst = sel_obj.inst()
  sel_cell = sel_inst.cell
  
  if len(sel_obj.path) == 1: #different behaviour when sel_obj.path length is 1 or is longer
    cplx_trans = sel_inst.trans #when length is 1, cell is in top level, we need to get its translation directly
  else:
    cplx_trans =  pya.ICplxTrans()
    for ie in sel_obj.path:
      #print(ie.specific_cplx_trans(), ie) #if the instance is part of some other cell, then all the translations are included in the ie_path
      cplx_trans = cplx_trans*ie.specific_cplx_trans() #accumulated transform

  inst_poly_dict = {}
  
  #check polygons for each layer inside the cell
  for idx in layout.layer_indices():
    
    region_base = pya.Region(sel_cell.begin_shapes_rec(idx), cplx_trans) #the base set of polygons for a particular layer
    
    #check if we are dealing with an instance array
    if not sel_inst.is_regular_array():
      if region_base.count()>0:
        #print(ly.get_info(idx).to_s(), region.count())
        #for poly in region.each_merged():
          #print(poly.to_s())
        inst_poly_dict[idx] = list(region_base.each_merged()) #if not an array, just consider the region on its own
    else: 
      #otherwise, consider the polygons over the regular array
      a, b = sel_inst.a, sel_inst.b
      na, nb = sel_inst.na, sel_inst.nb
      
      for ii in range(na):
        for jj in range(nb):
          
          region_ii_jj = pya.Region()
          region_ii_jj += region_base
          region_ii_jj.move(ii*a + jj*b)
          if region_ii_jj.count()>0:
        #print(ly.get_info(idx).to_s(), region.count())
        #for poly in region.each_merged():
          #print(poly.to_s())
            if ii == 0 and jj == 0: #initialize the list on the first index
              inst_poly_dict[idx] = list(region_ii_jj.each_merged())
            else: #add to the list subsequently
              inst_poly_dict[idx] += list(region_ii_jj.each_merged())
  return inst_poly_dict
  
def get_shape_poly_dict(sel_obj):
  '''
      Extracts polygons of a selected object known to be a shape into a <dict>

      Arguments
      --------

      :sel_obj: <pya.ObjectInstPath> Selected object

      Returns
      -------
      :poly_dict: <dict> of keys <int> layer index with values <list> of <pya.Polygon> for the associated layer
      '''
  lv = pya.Application.instance().main_window().current_view()

  sel_shape = sel_obj.shape
  cplx_trans =  pya.ICplxTrans()
  for ie in sel_obj.path:
    cplx_trans = cplx_trans*ie.specific_cplx_trans() #accumulated transform

  #works for box, polygon, path
  if is_pcell_guidingshape(sel_obj):
    layer_index = None
  else:
    layer_index = sel_shape.layer

  polygon = sel_shape.polygon
  polygon.transform(cplx_trans)
  poly_dict = {layer_index:[polygon]}
  return poly_dict


def sel_obj_abs_cplx_trans(sel_obj):
  '''
  Get the accumulated transform of a selected object <pya.ObjectInstPath>
  '''
  cplx_trans =  pya.ICplxTrans()
  for ie in sel_obj.path:
    cplx_trans = cplx_trans*ie.specific_cplx_trans() #accumulated transform
  return cplx_trans

def get_roundpath_poly_dict(sel_obj):
  '''
      Extracts polygons of a selected object known to be a ROUND_PATH into a <dict>

      Arguments
      --------
      kwargs
      :sel_obj: <pya.ObjectInstPath> Selected object

      Returns
      -------
      <dict> of keys <int> layer index with values <list> of <pya.Polygon> for the associated layer
      '''
  lv = pya.Application.instance().main_window().current_view()
  layout = lv.active_cellview().layout()

  #get polygon coordinates with respect to the current cellview
  cplx_trans =  pya.ICplxTrans()
  for ie in sel_obj.path:
    cplx_trans = cplx_trans*ie.specific_cplx_trans() #accumulated transform

  source_roundpath_inst = sel_obj.inst()
  trans = source_roundpath_inst.trans  
  pcell_parameters = source_roundpath_inst.pcell_parameters_by_name()
  
  dpath = pcell_parameters['path']
  #width = pcell_parameters['width']
  layer_info = pcell_parameters['layer']
  
  layer_index = layout.find_layer(layer_info)
  
  #dpath = selected_inst.pcell_parameters_by_name()['path']
  path = pya.Path()
  path.width = dpath.width/layout.dbu
  path.points = [pya.Point(pp.x/layout.dbu + trans.disp.x, pp.y/layout.dbu + trans.disp.y) for pp in dpath.each_point()] 
  polygon = path.polygon()
  polygon =  polygon.transform(cplx_trans)
  return {layer_index:[polygon]}
  
def _classify_selection(sel_obj):
  '''
  Classify a selected object with a string tag
  '''
  tag = None
  if sel_obj.is_cell_inst():
    tag = "instance"
    inst = sel_obj.inst()
    if inst.cell.name == "ROUND_PATH": #is a ROUND_PATH
      tag = "roundpath"
      
  else: #has shape
    tag = 'shape'
    if sel_obj.shape.is_path():
      tag = "path"
  return tag

def get_poly_dict(sel_obj, return_tag = False):
  '''
      Extracts polygons of a selected object into a <dict>

      Arguments
      --------

      :sel_obj: <pya.ObjectInstPath> Selected object

      kwargs
      :return_tag: <bool> Flag for returning the classification of the selected object as a string

      Returns
      -------
      :result: <dict> of keys <int> layer index with values <list> of <pya.Polygon> for the associated layer

      if :return_tag:
        :tag: <string> classification of :sel_obj:
      '''
  tag = _classify_selection(sel_obj)

  result = None
  if tag in ['shape', 'path']:
    result = get_shape_poly_dict(sel_obj)
  elif tag == 'roundpath': #just getting polygons as instance work better
    result = get_inst_poly_dict(sel_obj)#get_roundpath_polygon(sel_obj)
  elif tag == 'instance':
    result = get_inst_poly_dict(sel_obj)
   
  if return_tag:
    return result, tag
  else:
    return result

def valid_obj(sel_obj):
  '''Check that selected object is a cell inst or has a valid shape'''
  if sel_obj.is_cell_inst():
      return True
  else:
      #should have a shape
    if sel_obj.shape.is_valid():
        return True
  return False


def _valid_edge(edge):
  '''
  Filter for edges to be considered when decomposing from a polygon
  '''
  dbu = pya.CellView.active().layout().dbu
  edge_min_length = int(round(float(pya.Application.instance().get_config("snap-edge-min")) / dbu))
  return edge.length() >= edge_min_length


def poly_dict_to_edge_list(poly_dict, return_layer=True, filter_min_length=True, filter_manhattan=True, warn=True):
  '''

  Decomposes polygons from :poly_dict:, into edges, with options for filtering out particular edges

  (optional) Filter out edges below minimum length
  (optional) Filters out edges not on manhattan direction

  Arguments
  ---------
  :poly_dict: <dict> with layer index of the polygon as keys and <pya.Polygon> values
  :tag: <string> label of the originating object type

  kwargs
  :return_layer: <bool> Flag for returning the associated layers of the returned :edges:
  :filter_min_length: <bool> Filter out edges below a minimum length
  :filter_manhattan: <bool> Filter out edges not on a manhattan direction
  :warn: <bool> If no valid edges, emit a warning


  Returns
  -------

  :edges: <list> of <pya.Edges> decomposed from polygons from :poly_dict:

  if :return_layer:
    :layer_idxs: <list> of <int> layer indices corresponding to polygon the edge in :edges: belongs to

  '''
  if filter_min_length: #define variable for filtering out edge_min_length
    dbu = pya.CellView.active().layout().dbu
    edge_min_length = int(round(float(pya.Application.instance().get_config("snap-edge-min")) / dbu))

  # get visible layers
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  layer_prop_list = list(lv.each_layer())
  visible_layer_dict = dict([(layer_prop.layer_index(), layer_prop.visible) for layer_prop in layer_prop_list])

  # filter out all layers that are not currently visible
  poly_dict_new = {}
  for k in poly_dict.keys():
    if visible_layer_dict[k]:
      poly_dict_new[k] = poly_dict[k]
  poly_dict = poly_dict_new

  # check if any of the polygons in poly_dict are ports
  ports_indices = ports.layer_indices()
  inst_ports_indices = list(set(poly_dict.keys()) & set(ports_indices))
  ports_in_inst = (len(inst_ports_indices) > 0)  # any port indices in the poly dict?

  # if ports are present, process based only on the ports
  if ports_in_inst:
    port_polygons = []
    for idx in inst_ports_indices:
      port_polygons += poly_dict[idx]

    # use current selected layer as index
    sel_ly_info = ly.get_info(
      ly.find_layer(lv.current_layer.current().source_layer, lv.current_layer.current().source_datatype))
    sel_ly_idx = ly.layer(sel_ly_info.layer,
                          sel_ly_info.datatype)  # this defaults to the first available layer index when nothing is selected

    # resolve them to edges
    edges = []
    layer_idxs = []
    for pp in port_polygons:
      point, edge = ports.resolve_port(pp, False)

      conditions = [True]
      if filter_min_length:
        if edge.length() >= edge_min_length:
          conditions.append(True)
        else:
          conditions.append(False)
      if filter_manhattan:
        if edge.dx() * edge.dy() == 0:
          conditions.append(True)
        else:
          conditions.append(False)

      if all(conditions):
        edges.append(edge)
        layer_idxs.append(sel_ly_idx)

    if warn:
      if len(edges) == 0:
        pya.MessageBox.warning("Selection invalid",
                               "No valid edges in selected object.\nEdges too small or not manhattan.",
                               pya.MessageBox.Ok)
    if return_layer:
      return edges, layer_idxs
    else:
      return edges

  else:
    # if there are no ports, return all edges wider than snap-edge-min
    # create list of all the polygons

    edges = []
    layer_idxs = []

    for k in poly_dict:
      for poly in poly_dict[k]:  # go through each list in the dictionary
        for edge in poly.each_edge():

          conditions = [True]
          if filter_min_length:
            if edge.length() >= edge_min_length:
              conditions.append(True)
            else:
              conditions.append(False)
          if filter_manhattan:
            if edge.dx() * edge.dy() == 0:
              conditions.append(True)
            else:
              conditions.append(False)

          if all(conditions):
            edges.append(edge)
            layer_idxs.append(k)
    if warn:
      if len(edges) == 0:
        pya.MessageBox.warning("Selection invalid",
                               "No valid edges in selected object.\nEdges too small or not manhattan.",
                               pya.MessageBox.Ok)
    if return_layer:
      return edges, layer_idxs
    else:
      return edges


def get_edge_list(sel_obj, return_layer=True, filter_min_length=True, filter_manhattan=True):
  '''

  Gets polygons from :sel_obj: and decomposes them into edges, with options for filtering out particular edges

  (optional) Filter out edges below minimum length
  (optional) Filters out edges not on manhattan direction

  Arguments
  ---------
  :sel_obj: <pya.ObjectInstPath> Selected object

  kwargs
  :return_layer: <bool> Flag for returning the associated layers of the returned :edges:
  :filter_min_length: <bool> Filter out edges below a minimum length
  :filter_manhattan: <bool> Filter out edges not on a manhattan direction


  Returns
  -------

  :edges: <list> of <pya.Edges> decomposed from polygons from :poly_dict:

  if :return_layer:
    :layer_idxs: <list> of <int> layer indices corresponding to polygon the edge in :edges: belongs to

  '''
  poly_dict = get_poly_dict(sel_obj, return_tag = False)
  return poly_dict_to_edge_list(poly_dict, return_layer=return_layer, filter_min_length=filter_min_length, filter_manhattan=filter_manhattan)


def is_pcell_guidingshape(sel_obj):
  '''
  Checks if :sel_obj:, a <pya.ObjectInstPath>, is a PCell guiding shape
  '''
  if len(sel_obj.path)>0:
    sel_ie = sel_obj.path[-1]
    sel_inst = sel_ie.inst()
    if sel_inst.is_pcell() or sel_inst.cell.is_pcell_variant():
      return True
  return False

  #return None #do nothing if not a valid object




#####
# Math and Geometry
#####

def abs_min(values):
  '''
  Returns the smallest absolute value of a <list> of values
  '''
  #return the value closest to 0, regardless of sign
  min_dist = None
  min_val = None
  for v in values:
    if min_dist is None:
      min_dist = abs(v)
      min_val = v
    else:
      if abs(v)<min_dist:
        min_dist = abs(v)
        min_val = v
  return min_val
  

def sgn(val):
  '''
  sign implementation with case val = 0 returning 0

  Arguments
  ---------
  :val: <int> or <float> Value to obtain sign from

  Returns
  -------
  :v_s: -1, 0, or 1, sign of :val:
  '''
  v_s = 0
  if abs(val)>0:
    v_s = val/abs(val)
  return v_s


def project_point(e1, e2, p):
  '''
  Projects point :p: onto the line going through :e1: and :e2:

  Arguments
  ---------
  :e1: :e2: <pya.Point> Points defining the line to be projected on
  :p: <pya.Point> Point to be projected onto the line

  Return
  ------
  <pya.Point> point projected onto the line
  '''
  #edge point 1, edge point 2, point to project
  #assume all pya.Point

  v = e2 - e1#[e2.x - e1.x, e2.y - e1[1]]
  r = p - e1#[p[0] - e1[0], p1[1] - e1[1]]
  r_dot_v = v.sprod(r)#v[0]*r[0] + v[1]*r[1]
  v_abs = v.length()#norm(v) #(v[0]**2 + v[1]**2)**0.5
  t = r_dot_v/v_abs
  return e1 + v*(t/v_abs)#[v[0]*t/v_abs, v[1]*t/v_abs]




def intersect_lines(a0, a1, b0, b1):

  '''
  Intersection between two lines in 2D, y1(u) = a0 + a1*u and y2(v) = b0 + b1(v)

  a0 + a1*u = b0 + b1*v

  a0 - b0 = [a1 b1][-u v]T

  [a1 b1]^-1 (a0 - b0) = [-u v]T

  Arguments
  ---------
  :a0: <pya.Point>
  :a1: <pya.Point> or <pya.Vector>

  :b0: <pya.Point>
  :b1: <pya.Point> or <pya.Vector>

  Return
  ------
  <pya.Point> Intersection point

  '''
  
  a1x, a1y = a1.x, a1.y
  
  a0x = a0.x
  a0y = a0.y
  
  b0x = b0.x
  b0y = b0.y
  
  b1x = b1.x
  b1y = b1.y
  
  det = (a1x*b1y) - (a1y*b1x)
  cx = a0x - b0x
  cy = a0y - b0y
  u = -(b1y*cx - b1x*cy)/det
  
  return pya.Point(a0x+(a1x*u), a0y+(a1y*u))


def relative_coordinates(x, y, ie_path):
  '''
  Map (:x:,:y:) to relative coordinates corresponding to an instance within a hierarchy of instances,
  as specified by :ie_path:. Return the corresponding relative coordinates (:x_r:, :y_r:).

  Arguments
  ---------
  :x: :y: <int> <float> Coordinate values
  :ie_path: <list> of <pya.InstElement>, representing the hierarchy of instances of interest

  Return
  ------
  :x_r:, :y_r: <int> <float> Coordinate values

  '''
  x_r, y_r = x, y
  for ie in ie_path:
    x_r, y_r = coordinate_system_transform(x_r, y_r, ie.specific_cplx_trans())
  return x_r, y_r

def relative_coordinates_point(p, ie_path):
  '''
  Map :p: to relative coordinates corresponding to an instance within a hierarchy of instances,
  as specified by :ie_path:. Return the corresponding relative coordinates :p_r:.

    see relative_coordinates

  Too much work to set up multidispatch, and there is only 1 other use case anyways
  '''

  return type(p)(*relative_coordinates(p.x, p.y, ie_path))

def coordinate_system_transform(x, y, icplxtrans):
  '''
  Return coordinates corresponding to (:x:,:y:) if the coordinate axes were transformed according to transform :icplxtrans:

  Arguments
  ---------
  :x: :y: <int> <float> Origianl coordinate values
  :icplxtrans: <pya.ICplxTrans>, representing the new coordinate system

  Return
  ------
  :x_r:, :y_r: <int> <float> Coordinate values in the new coordinate system

  '''

  #Map current point into an instance's coordinate system
  #x, y: coordinate in base coordinate system
  #icplxtrans: new coordinate system, expressed as a transform on base coordinate system
  x0a = icplxtrans.disp.x
  y0a = icplxtrans.disp.y
  
  xp = x - x0a
  yp = y - y0a
  
  r = icplxtrans.mag
  th_rad = icplxtrans.angle*math.pi/180.

  if not icplxtrans.is_mirror():
    u_a =(r*math.cos(th_rad), r*math.sin(th_rad))
    v_a = (-r*math.sin(th_rad), r*math.cos(th_rad))
  else:
    u_a = (r*math.cos(th_rad), r*math.sin(th_rad))
    v_a = (r*math.sin(th_rad), -r*math.cos(th_rad))
    
  x_a = xp*u_a[0] + yp*u_a[1]
  y_a = xp*v_a[0] + yp*v_a[1]
  
  return x_a, y_a #coordinates in new coordinate system

def coordinate_system_transform_point(p, icplxtrans):
  '''
  Return coordinates p if the coordinate axes were transformed according to transform :icplxtrans:

  see coordinate_system_transform

  Too much work to use multidispatch, and there is only 1 other use case anyways
  '''
  return type(p)(*coordinate_system_transform(p.x, p.y, icplxtrans))

def rotate_axis(x, y, icplxtrans):
  '''
  Rotate coordinates corresponding to (:x:,:y:) according to transform :icplxtrans:

  Arguments
  ---------
  :x: :y: <int> <float> Original coordinate values
  :icplxtrans: <pya.ICplxTrans>, representing the new coordinate system

  Return
  ------
  :x_r:, :y_r: <int> <float> Rotated coordinates

  '''
  ang_rad = math.pi*icplxtrans.angle/180.
  
  sgn = 1
  if icplxtrans.is_mirror():
    sgn = -1
  
  r11 = math.cos(ang_rad)
  r12 = math.sin(ang_rad)
  r21 = -math.sin(ang_rad)
  r22 = math.cos(ang_rad)
  xr = round(r11*x + r12*y)
  yr = round(r21*x + r22*y)
  return xr, sgn*yr
  
#####
# Proximity Comparisons
#####

def closest_points(a_point_list, b_point_list):
  '''
  Brute-force closest pair of points from two sets of points :a_point_list: and :b_point_list:

  Arguments
  ------
  :a_point_list: <list> of <pya.Point>
  :b_point_list: <list> of <pya.Point>

  Returns
  ------
  :idx_a: <int> Index of point in :a_point_list: closest to points in :b_point_list:
  :idx_b: <int> Index of point in "b_point_list: closest to points in :a_point_list:
  :min_dist: <float> Minimum distance between the closest points of the two sets of points
  :min_dist_trans: <pya.Vector> The translation vector between the closest points

  '''
  min_dist = None
  min_dist_trans = None
  idx_a = 0
  idx_b = 0
  for i, a_point in enumerate(a_point_list):
    for j, b_point in enumerate(b_point_list):
      trans = (a_point - b_point)
      dist = trans.length()
      if min_dist is None or dist<min_dist:
        min_dist = dist
        min_dist_trans = trans
        idx_a = i
        idx_b = j

  return idx_a, idx_b, min_dist, min_dist_trans


def closest_path_segment_to_point(p_c, path):
  '''
    Identify the segment of a path closest to a point :p_c:

    Arguments
    ---------
    :p_c: <pya.Point> point of interest
    :path: <list> of <pya.Point> List of points constituting a path

    Return
    ------
    :idx_min: <int> index of :path: selected as the closest segment
    :point_min: <pya.Point> point projected on the path on the closest segment
    :end_flag: <bool> Flag determining if the selected segment is at the ends of the path
    :t_min: <float> projection parameterization (0<= t <= 1) on the path
    :val_min: <int> or <float> distance from :p_c: to the projected point on the closest path segment

    '''
  path_points = list(path.each_point())
  path_pairs = list(zip(path_points[:-1], path_points[1:]))

  idx_min = None
  val_min = None
  point_min = None
  t_min = None

  for i in range(len(path_pairs)):

    pp = path_pairs[i]
    p_start = pp[0]
    p_stop = pp[1]

    # calculate projection onto the segment
    v = p_stop - p_start
    v_c = p_c - p_start
    t = v.sprod(v_c) / (v.length() ** 2)

    # Check if the point lies on the current edge or not.
    # Else, restrict it to the endpoints of the segment
    if t > 1:
      t = 1
    elif t < 0:
      t = 0

    # projected point
    p_c_proj = p_start + t * v

    cur_length = (p_c_proj - p_c).length()

    if (val_min is None) or val_min > cur_length:
      idx_min = i
      val_min = cur_length
      point_min = p_c_proj
      t_min = t

  end_flag = False
  # check if at the endpoints of the path
  if idx_min == 0 and t_min == 0:
    end_flag = True
  elif idx_min == (len(path_pairs) - 1) and t_min == 1:
    end_flag = True

  return idx_min, point_min, end_flag, t_min, val_min


def closest_edges(a_edge_list, b_edge_list, filter_parallel=True, method = 'center'):
  '''
  Brute-force closest pair of edges from two sets of points :a_edge_list: and :b_edge_list:

  (Options)
  Considering only parallel edges

  Arguments
  ------
  :a_edge_list: <list> of <pya.Edge>
  :b_edge_list: <list> of <pya.Edge>

  :filter_parallel: <bool> Flag to consider only closest pairs of edges which are parallel
  :method: <string> center (default): compare center points, max: compare maximum distances between edges

  Returns
  ------
  :idx_a: <int> Index of point in :a_point_list: closest to points in :b_point_list:
  :idx_b: <int> Index of point in "b_point_list: closest to points in :a_point_list:
  :min_dist: <float> Minimum distance between the closest points of the two sets of points
  :min_dist_trans: <pya.Vector> The translation vector between the closest points

  '''

  min_dist = None
  min_dist_trans = None
  idx_a = 0
  idx_b = 0

  for i, a_edge in enumerate(a_edge_list):

    for j, b_edge in enumerate(b_edge_list):

      conditions = [True]
      if filter_parallel:
        va = a_edge.p2 - a_edge.p1
        vb = b_edge.p2 - b_edge.p1

        if va.length() * vb.length() == 0:  # if one of the edges is a point, continue and treat as point
          conditions.append(True)
        elif va.vprod(vb) == 0:  # is parallel
          conditions.append(True)
        else:
          conditions.append(False) #not parallel

      if all(conditions):


        if method.lower() == 'max':
          trans_list = [a_edge.p1 - b_edge.p1,
                        a_edge.p1 - b_edge.p2,
                        a_edge.p2 - b_edge.p1,
                        a_edge.p2 - b_edge.p2]
          dist_list = [vv.length() for vv in trans_list]

          idx_max = max(range(len(trans_list)), key=lambda x: dist_list[x])
          dist = dist_list[idx_max]
          trans = trans_list[idx_max]
        else: #default
          # center distance between points
          a_point = (a_edge.p1 + a_edge.p2) / 2.
          b_point = (b_edge.p1 + b_edge.p2) / 2.
          trans = (a_point - b_point)
          dist = trans.length()



        # if distance is smaller, update
        if min_dist is None or dist < min_dist:
          min_dist = dist
          min_dist_trans = trans
          idx_a = i
          idx_b = j

  return idx_a, idx_b, min_dist, min_dist_trans


def closest_edge_to_cursor(sel_obj, return_full=False):
  '''
    Closest edge of selected object :sel_obj: with cursor
    (Options) Returns also cursor location and layer info of the edge
.
    Arguments
    ------
    :sel_obj: <pya.ObjectInstPath> Selected object

    kwargs
    :return_full: <bool> Flag to return also :cursor_point: and :layer_idx:

    Returns
    ------
    :fixed_edge: <pya.Edge> Edge from the selected object closest to the curent cursor position

    if return_full:
        :cursor_point: <pya.Point> Current cursor position (in dbu, relative to the current cellview)
        :fixed_layer: <int> Layer index of polygon the edge originates from
    '''
  lv = pya.Application.instance().main_window().current_view()

  # Obtain the relative position with respect to the current cellview of the cursor, not snapped to grid point
  cursor_point = cursor_position(in_dbu=True, relative=True, on_grid=False)

  # get a list of the points of the polygon of the shape
  # assume fixed_obj is a box, polygon, or path. Will convert to polygon
  # fixed_points, fixed_edges = cursor_edgeenter_points(fixed_obj)

  # get the edges of the selected object
  fixed_edges, fixed_layer_idx = get_edge_list(sel_obj, return_layer=True, filter_min_length=True, filter_manhattan=True)
  if len(fixed_edges)==0:
    return None
  fixed_edges_center_points = [(ee.p1 + ee.p2)/2. for ee in fixed_edges]
  idx_fixed, idx_move, min_dist, min_dist_trans = closest_points(fixed_edges_center_points, [cursor_point])

  fixed_edge = fixed_edges[idx_fixed]
  fixed_layer = fixed_layer_idx[idx_fixed]
  if return_full:
    return fixed_edge, cursor_point, fixed_layer
  return fixed_edge  # return only the edge

#####
# Path Manipulation
#####


def segment_path_corners(path, radius):
  '''
  Segment a manhattan path with corners of certain radius

  Arguments
  ---------
  :path: <list> of <pya.Point> a Manhattan path
  :radius: <int> <float> Radius of the corners to be segmented (depending on whether path is in dbu or in measurement units)

  Return
  ------
  :segments: <list> of <list> of <pya.Point>, path segments between corners

  :corners: <list> of <list> of <pya.Point>, each corner is  a list [p0, p1, p2], p1 is the
  corner, p0 and p2 are ends of the corner bend

  :boxes: <list> of <list> of <pya.Point>, corners that do not satisfy clearance requirements. each corner is  a
  list [p0, p1, p2], p1 is the corner, p0 and p2 are ends of the corner bend

  '''
  # points of the current path
  points = [pp for pp in path.each_point()]

  pp_prev = None
  pp_cur = None
  pp_next = None

  segments = []
  seg_cur = []

  corners = []
  boxes = []

  # crawl along the path
  for i in range(len(points)):
    # print('----')
    # print(pp_prev)
    # print(pp_cur)
    # print(pp_next)

    pp_cur = points[i]
    if i + 1 < len(points):
      pp_next = points[i + 1]
    else:
      pp_next = None

    # check if this is a corner
    if not (pp_prev is None or pp_next is None):
      v_p = pp_cur - pp_prev
      v_n = pp_next - pp_cur
      is_corner = ((abs(v_p.sprod(v_n)) == 0))
      has_space = (v_p.length() >= radius) and (v_n.length() >= radius)  # path and radius should be in same units
    else:
      is_corner = False

    pp_prev = pp_cur
    # print(is_corner)

    if is_corner:
      if has_space:
        # if it's a corner segment, draw a polygon indicating the corner when possible
        # this is a valid corner
        p1 = pp_cur
        p0 = p1 - radius * (v_p / v_p.length())
        p2 = p1 + radius * (v_n / v_n.length())

        # poly = pya.Polygon([p0, p1, p2])
        # topcell.shapes(layer_idx).insert(poly)
        corners.append([p0, p1, p2])  # p1 is the corner

        # adjust the next previous value
        pp_prev = p2
        # make a new segment
        seg_cur.append(p0)
        segments.append(seg_cur)
        seg_cur = [p2]
      else:

        p1 = pp_cur
        p0 = p1 - type(p1)(radius, radius)
        p2 = p1 + type(p1)(radius, radius)

        # box = pya.Box(p0, p2)#
        # topcell.shapes(layer_idx).insert(box) #shows the clearance for the path
        boxes.append([p0, p1, p2])
        # don't make a new segment
        seg_cur.append(p1)
    else:
      seg_cur.append(pp_cur)

  segments.append(seg_cur)  # append the final segment

  return segments, corners, boxes


def simplify_path(path_points):
  '''
  Remove redundant points on the path
  '''
  result = [path_points[0]]
  for i in range(1, len(path_points)-1):
    if not( (path_points[i-1].x == path_points[i].x == path_points[i+1].x) or (path_points[i-1].y == path_points[i].y == path_points[i+1].y) ):
      result.append(path_points[i])
  result.append(path_points[-1])
  return result
