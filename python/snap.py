'''
snap.py

snapper v 0.0.1
Copyright 2023 Jason Mak

This file is part of snapper.

snapper is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

snapper is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with snapper. If not, see <https://www.gnu.org/licenses/>.

'''

import pya

import snapper.utils as utils
import snapper.pcell as pcell


def snap_obj(orientation = 0):
  '''
Snaps a selected object to 1) Cursor 2) an object 3) aligned to multiple objects

Align by :orientation: (0: center, -1: to left/bottom, 1: to top/right)



Assumptions
----------
All snapping occurs on a Manhattan basis. All edges that are snapped are vertically or horizontally aligned. Behaviour of rotated edges not guaranteed.

If selected object is an instance, we assume it either has :ports: drawn at its ports, or we will assume we want to snap only to edges on the boundaries of its bounding box.


Arugments
---------
:orientation: Sets the alignment of the snapping edge (0: center, -1: to left/bottom, 1: to top/right)


GUI Inputs
----------

Situation 1: 1 selected object
Situation 2: 2 selected objects
Situation 3: 3 or more selected objects.


Behaviour
----------
Situation 1:
The edge of the selected object closest to the cursor will be moved to the grid point closest to the cursor by translation (no rotation or mirroring)

Situation 2:
The first selected object is brought to the second object, to join the closest edges.

Situation 3:
The first selected object is brought to the second object, for the axis traverse to its edge, and then the same to the third selected object.
e.g. For aligning a 90 degree bend to 2 waveguides in orthogonal directions.

  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell
  
  
  object_selection = list(lv.object_selection) #important to ensure ordering stays as expected, or else may behave like a set?
  object_selection = sorted(object_selection, key =lambda x: x.seq)

  ### CASE: 1 selected objects. Move the first selection to the cursor location
  if len(object_selection) == 1:
    #### Assess current state of cursor position and 
    
    #cursor position, in dbu
    cursor_point = utils.cursor_position(in_dbu=True, relative=True, on_grid=True)
  
    move_obj = object_selection[0]
    move_edges = utils.get_edge_list(move_obj, return_layer=False, filter_min_length=True, filter_manhattan=True)
    if len(move_edges) == 0:
      return
    #select the edge with its center point closest to the cursor
    move_edge_centers = [(ee.p1 + ee.p2)/2. for ee in move_edges]
    idx_move, idx_fixed, min_dist, center_trans = utils.closest_points(move_edge_centers, [cursor_point])


    #based on the orientation, construct the align vector
    move_edge = move_edges[idx_move]
    center_trans = cursor_point - move_edge_centers[idx_move]

    move_edge_width = move_edge.length()
    align_shift = -move_edge_width/2.
    
    if move_edge.dx() == 0:
      align_vector = pya.Vector(0, align_shift*orientation)
    elif move_edge.dy() == 0:
      align_vector = pya.Vector(align_shift*orientation, 0)
    else:
      align_vector = pya.Vector(0,0)

    #check if move_obj is inside an instance
    trans_vector = center_trans + align_vector
    #when using snap-align mode, and the translation is two axes, translate along edge axis first
    #if (str(pya.Application.instance().get_config('snap-align')).lower() == "true") and not trans_vector.x*trans_vector.y == 0: #obsolete?
    #  move_edge_norm_vector = move_edge.d()/move_edge.length()
    #translate only along direction of edge
    #  trans_vector = pya.Vector(trans_vector.x*abs(move_edge_norm_vector.x), trans_vector.y*abs(move_edge_norm_vector.y))
    trans = pya.Trans(trans_vector)  
  
    lv.transaction("snap obj - 1")
    
    if len(move_obj.path)>0:
      #if inside an instance, move the highest level instance
      move_obj.path[0].inst().transform(trans)
      lv.object_selection = []
      
    else:
      move_obj.shape.transform(trans) 
      
        
    lv.commit()
    
    
  ### CASE: 2 selected objects. Move the first selection to the second selection
  elif len(object_selection) == 2:
    
    #The most recent selection is added to the beginning of the list
    fixed_obj = object_selection[1] #last selected object is first in the list to be fixed
    move_obj = object_selection[0] 
    
    #shapes: polygons, boxes, path: turn it into polygons, get all its edges
    #instances: check if it has ports. otherwise, use edges that intersect its bounding box (typically waveguides)
    
    #get the edges for fixed and move obj
    fixed_edges = utils.get_edge_list(fixed_obj, return_layer=False, filter_min_length=True,
                                        filter_manhattan=True)
    if len(fixed_edges) == 0:
      return
    move_edges = utils.get_edge_list(move_obj, return_layer=False, filter_min_length=True,
                                        filter_manhattan=True)
    if len(move_edges) == 0:
      return
    idx_fixed, idx_move, min_dist, min_dist_trans = utils.closest_edges(fixed_edges, move_edges, filter_parallel=True, method='center')

    #based on the orientation, construct the align vector
    fixed_edge = fixed_edges[idx_fixed]
    move_edge = move_edges[idx_move]
    center_trans = (fixed_edge.p1 + fixed_edge.p2)/2. - (move_edge.p1 + move_edge.p2)/2.
    
    fixed_edge_width = fixed_edge.length()
    move_edge_width = move_edge.length()
    align_shift = (fixed_edge_width - move_edge_width)/2.
    
    if move_edge.dx() == 0:
      align_vector = pya.Vector(0, align_shift*orientation)
    elif move_edge.dy() == 0:
      align_vector = pya.Vector(align_shift*orientation, 0)
    else:
      align_vector = pya.Vector(0,0)
    
    #trans = pya.Trans(min_dist_trans + align_vector)
    #check if move_obj is inside an instance
    trans_vector = center_trans + align_vector
    
    #when using snap-align mode, and the translation is two axes, translate along edge axis first
    #if (str(pya.Application.instance().get_config('snap-align')).lower() == "true") and not trans_vector.x*trans_vector.y == 0:
    #  move_edge_norm_vector = move_edge.d()/move_edge.length()
    #translate only along direction of edge
    #  trans_vector = pya.Vector(trans_vector.x*abs(move_edge_norm_vector.x), trans_vector.y*abs(move_edge_norm_vector.y))
    
    trans = pya.Trans(trans_vector)    
    
    ###### Transaction
    lv.transaction("snap obj - 2")
    
    if len(move_obj.path)>0:
      #if inside an instance, move the highest level instance
      move_obj.path[0].inst().transform(trans)
      lv.object_selection = []
      
    else:
      move_obj.shape.transform(trans) 
      
        
    lv.commit()
        
  
  ### CASE:  multiple selected objects. Default in -align mode. Partially align successively to each selected fixed object.
  elif len(object_selection) > 2:
    
    #The most recent selection is added to the beginning of the list
    fixed_obj_list = object_selection[1:] #last selected object is first in the list to be fixed
    move_obj = object_selection[0] 
    
    #shapes: polygons, boxes, path: turn it into polygons, get all its edges
    #instances: check if it has ports. otherwise, use edges that intersect its bounding box (typically waveguides)

    move_edges = utils.get_edge_list(move_obj, return_layer=False, filter_min_length=True,
                                        filter_manhattan=True)
    if len(move_edges) == 0:
      return
    ###### Transaction
    lv.transaction("snap obj - multi")
    
    for fixed_obj in fixed_obj_list:
      #get the polygons _dict for the fixed object
      fixed_edges = utils.get_edge_list(fixed_obj, return_layer=False, filter_min_length=True,
                                           filter_manhattan=True)
      if len(fixed_edges) == 0:
        return

      #find closest point: compute all distance pairs between edges and find the minimum
      idx_fixed, idx_move, min_dist, min_dist_trans = utils.closest_edges(fixed_edges, move_edges, filter_parallel=True, method='center')
      #based on the orientation, construct the align vector
      fixed_edge = fixed_edges[idx_fixed]
      move_edge = move_edges[idx_move]
      center_trans = (fixed_edge.p1 + fixed_edge.p2) / 2. - (move_edge.p1 + move_edge.p2) / 2.
      
      fixed_edge_width = fixed_edge.length()
      move_edge_width = move_edge.length()
      align_shift = (fixed_edge_width - move_edge_width)/2.
      
      if move_edge.dx() == 0:
        align_vector = pya.Vector(0, align_shift*orientation)
      elif move_edge.dy() == 0:
        align_vector = pya.Vector(align_shift*orientation, 0)
      else:
        align_vector = pya.Vector(0,0)
      
      #trans = pya.Trans(min_dist_trans + align_vector)
      #check if move_obj is inside an instance
      trans_vector = center_trans + align_vector
      #when using -align mode, and the translation is two axes, translate along edge axis first
      #if True:#(str(pya.Application.instance().get_config('-align')).lower() == "true") and not trans_vector.x*trans_vector.y == 0:
      move_edge_norm_vector = move_edge.d()/move_edge.length()
      #translate only along direction of edge
      trans_vector = pya.Vector(trans_vector.x*abs(move_edge_norm_vector.x), trans_vector.y*abs(move_edge_norm_vector.y))
      trans = pya.Trans(trans_vector)    
    
  
    
      if len(move_obj.path)>0:
        #if inside an instance, move the highest level instance
        move_obj.path[0].inst().transform(trans)
        lv.object_selection = []
        
      else:
        move_obj.shape.transform(trans) 
      
        
    lv.commit()


#############################
# SNAP PATH
#############################

def __path_check_move_obj(move_obj):
  '''
  Checks the object intended to be moved to be a valid path/path related object, assigns a tag

  Arguments
  ---------
  :move_obj: <pya.ObjectInstPath> Selected object to be checked

  Returns
  -------
  :move_obj_type: <string> Tag, identifying the originating type of path
  :move_path: <pya.Path> Extracted <pya.Path> from :move_obj:

  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()

  move_obj_type = None

  if move_obj.is_cell_inst():  # make function compatible when selecting PCell by instances

    if pcell.is_path_pcell_inst(move_obj):  # inst_move.is_pcell(): #is a PCell, with a path
      move_obj_type = "PathPCell"
      # path along its internal coordinates
      move_path = pya.Path(pcell.get_path_pcell_path(move_obj) * round(1. / ly.dbu))  # in dbu

  elif move_obj.shape.is_valid():  # is not an instance, and is valid shape

    if utils.is_pcell_guidingshape(move_obj):  # is a PCell GuidingShape
      move_obj_type = "GuidePath"
      move_path = move_obj.shape.path

    elif move_obj.shape.is_path():  # is a Path
      move_obj_type = "Path"
      move_path = move_obj.shape.path

  if move_obj_type is None: #move_path will have been declared otherwise
    pya.MessageBox.warning("Selection not a path",
                           "First selection must be a path or supported PCell", pya.MessageBox.Ok)

    return None, None

  return move_obj_type, move_path


def _snap_path_end_point_move(path_points, idx_move, cursor_point, align_adjust, orientation):
  '''
  Adjust an end point of :path_points: (as determined by :idx_move:) to :cursor_point:,
  making the appropriate alignment based on :align_adjust: and :orientation:.

  Assumes the end point terminate on a manhattan direction (throws warning msgbox if not manhattan)

  Arguments
  ---------
  :path_points: <list> of <pya.Point> Points of a path
  :idx_move: <int> index of list :path_points: for the selected endpoint
  :cursor_point: <pya.Point> location that we want the path end point to relocate to
  :align_adjust: <float> adjustment factor to align path center with target center (could be an edge)
  :orientation: <float> alignment flag for determining left/center/right justification of alignment

  Returns
  -------
  :path_points_moved: <list> of <pya.Point> Adjusted points of the path
  '''
  if len(path_points) > 2:  # non-trivial line

    # determine the orientation of the path and the edge
    if idx_move == 0:
      # at the start of the path
      end_point = path_points[0]
      prev_point = path_points[1]
      prev_dir = path_points[2] - path_points[1]


    elif idx_move == len(path_points) - 1:  # move_path.points - 1: #at end of the path
      end_point = path_points[idx_move]
      prev_point = path_points[idx_move - 1]
      prev_dir = path_points[idx_move - 2] - path_points[idx_move - 1]

    # determine which direction corresponds to the translation
    path_edge = end_point - prev_point

    # assume manhattan connections for simplicity
    # assert(path_edge.x*path_edge.y == 0) #Manhattan
    if not path_edge.x * path_edge.y == 0:
      pya.MessageBox.warning("Path is not Manhattan",
                             "Snapping ends of path needs to be manhattan", pya.MessageBox.Ok)

      return None

    # move the end point and its previous point in an angle preserving way (along its current direction)
    if path_edge.x == 0:  # vertical
      align_adjust = pya.Vector(orientation * align_adjust, 0)
      end_point_dir = pya.Point(0, 1)

    elif path_edge.y == 0:  # horizontal
      # displace last point by x, second to last point by y
      align_adjust = pya.Vector(0, orientation * align_adjust)
      end_point_dir = pya.Point(1, 0)

    end_point_moved = cursor_point + align_adjust
    prev_point_moved = utils.intersect_lines(end_point_moved, end_point_dir, prev_point, prev_dir)

  else:  # for trivial line segment (2 points)
    # keep the line direction the same, but project the previous point onto the new position
    if idx_move == 0:
      # at the start of the path
      end_point = path_points[0]
      prev_point = path_points[1]
    elif idx_move == len(path_points) - 1:  # move_path.points - 1: #at end of the path
      end_point = path_points[idx_move]
      prev_point = path_points[idx_move - 1]

    path_edge = end_point - prev_point

    if not path_edge.x * path_edge.y == 0:
      pya.MessageBox.warning("Path is not Manhattan",
                             "Snapping ends of path needs to be manhattan", pya.MessageBox.Ok)

      return None

    # move the end point and its previous point in an angle preserving way (along its current direction)
    if path_edge.x == 0:  # vertical
      align_adjust = pya.Vector(orientation * align_adjust, 0)
      end_point_dir = pya.Point(0, 1)

    elif path_edge.y == 0:  # horizontal
      # displace last point by x, second to last point by y
      align_adjust = pya.Vector(0, orientation * align_adjust)
      end_point_dir = pya.Point(1, 0)

    end_point_moved = cursor_point + align_adjust
    prev_point_moved = utils.project_point(end_point_moved, end_point_moved + end_point_dir, prev_point)

  path_points_moved = path_points.copy()
  if idx_move == 0:
    # at the start of the path
    path_points_moved[0] = end_point_moved
    path_points_moved[1] = prev_point_moved
  elif idx_move == len(path_points_moved) - 1:  # at end of the path
    path_points_moved[idx_move] = end_point_moved
    path_points_moved[idx_move - 1] = prev_point_moved

  return path_points_moved  # end_point_moved, prev_point_moved

def snap_path(orientation = 0):
  '''

Snaps a selected path/PCell with path guiding shape/guiding shape path of a PCell
to 1) cursor a) at the end of the path, b) along the path 2) the end of the path to the [alignment] of an edge of a subsequently selected object(s)

:orientation: (0: center, -1: to left/bottom, 1: to top/right)


For simplicity, we refer to any path/PCell with path guiding shape/guiding shape path of a PCell generally as path

Arguments
--------
:orientation: (0: center, -1: to left/bottom, 1: to top/right)


Assumptions
----------
All snapping occurs on a Manhattan basis. All edges that are snapped are vertically or horizontally aligned. Behaviour of rotated edges not guaranteed.

If selected object is an instance, we assume it either has :ports: drawn at its ports, or we will assume we want to snap only to edges on the boundaries of its bounding box.

GUI Inputs
----------

Situation 1: 1 selected object (the path)
  a) Cursor close to the ends
  b) Cursor close to a segment along the path
Situation 2: 2 selected or more selected objects



Behaviour
----------
Situation 1:
a) The end of the path is extended/contracted such that the end of the path ends on the grid point closest to the cursor. Angles of the previous line segments are preserved in the translation of the end point.
b) The line segment of the path closest to the cursor is adjusted such that the segment now overlaps the grid point closest to the cursor. Angles of surrounding line segments are preserved.

Situation 2:
The selected path (first selection) has its end points translated to meet the edges of the subsequent selections.

See scripts
----------
snap_path
snap_path_align_bot_left
snap_path_align_top_right

  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell

  object_selection = list(lv.object_selection)  # important to ensure ordering stays as expected
  object_selection = sorted(object_selection, key=lambda x: x.seq)

  if len(object_selection) == 1:
    # Align path to cursor, if the path is the sole selection. Behaviour is cursor dependent.

    # get and process obj_mov
    move_obj = object_selection[0]  # a path, PCell with guiding path, or guiding path of a PCell

    move_obj_type, move_path = __path_check_move_obj(move_obj)  # classify and extract the relevant components of move_obj

    # cursor position, at nearest grid point, relative to the current coordinate system of cellview
    cursor_point = utils.cursor_position(in_dbu=True, relative=True, on_grid=True)

    # if shape is within any instance, use the coordinates of the same level as the path
    # This automatically accounts for guiding shapes
    if len(move_obj.path) > 0:
      cursor_point = utils.relative_coordinates_point(cursor_point, move_obj.path)
    if move_obj_type is "PathPCell":
      # if selecting a PCell, additionally translate the cursor into the PCell's coordinate system
      cursor_point = utils.coordinate_system_transform_point(cursor_point, move_obj.inst().cplx_trans)

    if move_obj_type is None:
      return None  #if not a recognized move_obj type, do nothing

    # get location on the path
    idx_min, point_min, end_flag, t_min, val_min = utils.closest_path_segment_to_point(cursor_point, move_path)

    if not end_flag:  # if in the middle of the path,  the path to the grid point closest to the cursor

      move_path_points = list(move_path.each_point()) #to list of points


      ## Iterate over path to get current segment
      prev_point = move_path_points[idx_min]
      next_point = move_path_points[idx_min + 1]
      current_direction = next_point - prev_point
      orthogonal_direction = pya.Point(current_direction.y, -current_direction.x)

      if len(move_path_points) > 2:

        # try to get prev2 and next2 points, but handle edge cases
        if idx_min == 0:
          # get directions
          next_direction = move_path_points[idx_min + 2] - next_point
          prev_direction = orthogonal_direction  # no previous point
        elif idx_min + 2 >= len(move_path_points):
          prev_direction = prev_point - move_path_points[idx_min - 1]
          next_direction = orthogonal_direction  # orthogonal to current direction, no next point
        else:
          # get directions
          prev_direction = prev_point - move_path_points[idx_min - 1]
          next_direction = move_path_points[idx_min + 2] - next_point

        # new points through intersection, extending current directions
        prev_point_new = utils.intersect_lines(cursor_point, current_direction, prev_point, prev_direction)
        next_point_new = utils.intersect_lines(cursor_point, current_direction, next_point, next_direction)

      else:  # 2 point line, simply move the line to intersect cursor
        correction = cursor_point - point_min
        prev_point_new = prev_point + correction
        next_point_new = next_point + correction

      move_path_points[idx_min] = prev_point_new
      move_path_points[idx_min + 1] = next_point_new

      new_path_points = move_path_points

      lv.transaction("snap path - path to grid")

      # handle the various kinds of supported path
      if move_obj_type is "PathPCell":
        new_path = pya.DPath(new_path_points, move_path.width)
        new_path = new_path * ly.dbu #adjust back to microns
        pcell.set_path_pcell_path(move_obj, new_path)  # pcell is pruned within this function call

        object_selection = [move_obj]  # create a new selection of the selected objects
        lv.object_selection = object_selection

      elif move_obj_type in ["Path", "GuidePath"]:

        # update move_obj.shape with the path with updated points
        move_obj.shape.path = pya.Path(new_path_points, move_path.width)
        cell.refresh()
        pcell.refresh_pcell_path(move_obj)  # update if is pcell guiding shape
      lv.commit()

    else:  # the end of the path to cursor

      # Cursor is in coordinates on the same frame of reference as the path
      # Hence, we can compare without taking the absolute coordinates of the path
      move_path_points = list(move_path.each_point())
      move_path_end_points = [move_path_points[0], move_path_points[-1]]
      move_path_end_points_idx = [0, len(move_path_points) - 1]

      # find closest point: compute all distance pairs and find the minimum
      idx_move_partial, idx_fixed, min_dist, min_dist_trans = utils.closest_points(move_path_end_points, [cursor_point])
      idx_move = move_path_end_points_idx[idx_move_partial]

      # we now know the surfaces closest to each other
      # we want to  to the right side


      align_adjust = - move_path.width / 2

      # determine end point and previous endpoint to move
      new_path_points = _snap_path_end_point_move(move_path_points, idx_move, cursor_point, align_adjust, orientation)
      if new_path_points is None:
        return  #exit when not valid

      lv.transaction("snap path - to other")

      # handle the various kinds of supported path
      if move_obj_type is "PathPCell":

        new_path = pya.DPath(new_path_points, move_path.width)
        new_path = new_path * ly.dbu
        pcell.set_path_pcell_path(move_obj, new_path)

        object_selection = [move_obj]  # create a new selection of the selected objects
        lv.object_selection = object_selection

      elif move_obj_type in ["Path", "GuidePath"]:

        # update move_obj.shape with the path with updated points
        move_obj.shape.path = pya.Path(new_path_points, move_path.width)
        cell.refresh()
        pcell.refresh_pcell_path(move_obj)  # update if is pcell guiding shape
      lv.commit()


  elif len(object_selection) >= 2:  # Selection of at least 2 objects
    move_obj = object_selection[
      0]  # First object: Path, PCell with path guiding shape, or a path guiding shape
    fixed_obj_list = object_selection[1:]  # a polygon

    ### CHECK move_obj type
    move_obj_type, move_path = __path_check_move_obj(move_obj)
    if move_obj_type is None:
      return #exit if incorrect object type

    move_path_abs = move_path.transformed_cplx(utils.sel_obj_abs_cplx_trans(move_obj))  # absolute path coordinates
    move_path_points = list(move_path.each_point())

    move_path_abs_points = list(
      move_path_abs.each_point())  # to make the comparison between the path and the fixed selections in accumulated transform coordinates
    move_path_abs_end_points = [move_path_abs_points[0], move_path_abs_points[-1]]
    move_path_abs_end_points_idx = [0, len(move_path_abs_points) - 1]

    for fixed_obj in fixed_obj_list:
      # get a list of the points of the polygon of the shape, in absolute coordinates
      fixed_edges, fixed_layer_idx = utils.get_edge_list(fixed_obj, return_layer=True, filter_min_length=True, filter_manhattan=True)

      if len(fixed_edges) == 0:
        return

      fixed_edges_centers = [(ee.p1 + ee.p2)/2. for ee in fixed_edges]
      # find closest point: compute all distance pairs and find the minimum

      idx_fixed, idx_move_partial, min_dist, min_dist_trans = utils.closest_points(fixed_edges_centers, move_path_abs_end_points)
      idx_move = move_path_abs_end_points_idx[idx_move_partial]

      # we now know the surfaces closest to each other
      # we want to  to the right side

      edge_width = fixed_edges[idx_fixed].length()
      align_adjust = (edge_width - move_path.width) / 2

      cursor_point = fixed_edges_centers[idx_fixed]  # set the selected closest point as the cursor location

      if len(move_obj.path) > 0:
        cursor_point = utils.relative_coordinates_point(cursor_point, move_obj.path)

      if move_obj_type is "PathPCell":
        # if selecting a PCell, additionally translate the cursor into the PCell's coordinate system
        cursor_point = utils.coordinate_system_transform_point(cursor_point,move_obj.inst().cplx_trans)

      # determine end point and previous endpoint to move
      new_path_points = _snap_path_end_point_move(move_path_points, idx_move, cursor_point, align_adjust, orientation)
      if new_path_points is None:
        pass  # return #exit when not valid

      # update points list
      move_path_points = new_path_points

      lv.transaction("snap path - to other")

      # handle the varioud kinds of supported path
      if move_obj_type is "PathPCell":

        new_path = pya.DPath(new_path_points, move_path.width)
        new_path = new_path * ly.dbu
        pcell.set_path_pcell_path(move_obj, new_path)

        object_selection = [move_obj]  # create a new selection of the selected objects
        lv.object_selection = object_selection

      elif move_obj_type in ["Path", "GuidePath"]:

        # update move_obj.shape with the path with updated points
        move_obj.shape.path = pya.Path(new_path_points, move_path.width)
        cell.refresh()
        pcell.refresh_pcell_path(move_obj)  # update if is pcell guiding shape

      lv.commit()

def snap_box(axis, direction):
  '''
  snap_box

Adjust an edge of the currently selected box to the nearest edge of the second selection,
based on :axis: and :direction:

Arguments
---------
:axis: <string> 'x' or 'y', in coordinate system of the box
:direction: <float> -1 or 1, along coordinate system of the box

Note: coordinate system of the box may not be the same as the current view coordinate system
Hence, need to use utils.orientation to convert the current view coordinates to the reference frame of the box

GUI Inputs
----------
2 selected objects:

1) a Box

2) an object

Behaviour
----------

The :axis:/:direction: (e.g. x, -1 == BOTTOM, y, +1 == RIGHT) edge of the currently selected box is
adjusted to overlap with the edge of the object closest to the box.

If the object is an instance, we consider its bounding box.


  '''
  lv = pya.Application.instance().main_window().current_view()
  layout = lv.active_cellview().layout()
  cell = lv.active_cellview().cell
  
  #grid, in dbu
  g_dbu = utils.grid()
  
  #cursor position
  p_c = utils.cursor_position(in_dbu=True, relative=True, on_grid=True)
  x_c_g, y_c_g = p_c.x, p_c.y
  #current selected objects
  #object_selection = lv.object_selection #this is an attribute
  object_selection = list(lv.object_selection) #important to ensure ordering stays as expected, or else may behave like a set?
  object_selection = sorted(object_selection, key=lambda x: x.seq)
  
  if len(object_selection)<1:
    pya.MessageBox.warning("Insufficient selection",
                                 "At least 1 object must be selected", pya.MessageBox.Ok)
    return  
  
  move_obj = object_selection[0] #check that this is a box
#  print("Is box:", move_obj.shape.is_box())
  if move_obj.is_cell_inst() or not move_obj.shape.is_box():
    pya.MessageBox.warning("Box not selected",
                                 "First selected object mus be box", pya.MessageBox.Ok)
    return  
  
  #compile list of reference points to  towards
  u_list = []
  if len(object_selection)>1:
    fixed_obj = object_selection[1] 
    fix_poly_dict = utils.get_poly_dict(fixed_obj) #get the polygon associated with the fixed selected object
    
    #based on the axis, simplify the polyon into list of values  
    for k, poly_list in fix_poly_dict.items():
      for poly in poly_list:
        for p in poly.each_point_hull(): #assume it's convex
          if axis == 'x':
            u_list.append(p.x)
          else:
            u_list.append(p.y)
  else:
    if axis == 'x':
      u_list.append(x_c_g)
    else:
      u_list.append(y_c_g)

  
  #modify the box
  box = move_obj.shape.box
  
  if axis == 'x' and direction == -1:      
    diff =  utils.abs_min([move_obj.shape.box.p1.x - u for u in u_list])
    box.p1 = pya.Point(box.p1.x - diff, box.p1.y)
   
  elif axis == 'x' and direction == 1:
    diff =  utils.abs_min([u - move_obj.shape.box.p2.x for u in u_list])
    box.p2 = pya.Point(box.p2.x + diff, box.p2.y)
    
  elif axis == 'y' and direction == -1:
    diff =  utils.abs_min([move_obj.shape.box.p1.y - u for u in u_list])
    box.p1 = pya.Point(box.p1.x, box.p1.y - diff)
   
  elif axis == 'y' and direction == 1:
    diff =  utils.abs_min([u - move_obj.shape.box.p2.y for u in u_list])
    box.p2 = pya.Point(box.p2.x, box.p2.y + diff)

  lv.transaction("snap box")
  move_obj.shape.box = box
  cell.refresh()
  lv.commit()

  return
    