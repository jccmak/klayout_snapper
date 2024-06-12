'''
pcell.py

snapper v 0.0.1
Copyright 2023 Jason Mak

This file is part of snapper.

snapper is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

snapper is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with snapper. If not, see <https://www.gnu.org/licenses/>.
'''

import pya
import math


pcell_names = ['ArrayPath', 'ROUND_PATH', 'ClearancePath']


def is_path_pcell_inst(sel_obj):
  '''
  Determine if :sel_obj: <pya.ObjectInstPath> is a PCell with a <pya.DPath> or <pya.Path> guiding shape
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  
  if sel_obj.is_cell_inst():
    #current selection is an instance
    sel_inst = sel_obj.inst()
    if sel_inst.is_pcell():
      #current selection is a pcell
      
      #check if its base cell is one of the known path based PCells
      if sel_inst.cell.name in pcell_names:
        return True
      else: #general handling
        pcell_parameters = sel_inst.pcell_parameters_by_name()
        for k,v in pcell_parameters.items():
          if type(v) is pya.Path:
            return True
          elif type(v) is pya.DPath:
            return True
  return False
  
def get_path_pcell_path(sel_obj):
  '''
  Return the <pya.Path> or <pya.DPath> guiding shape from :sel_obj: <pya.ObjectInstPath>
  '''
  #assume sel_obj is known to be a pcell
  sel_inst = sel_obj.inst()
  
  pcell_parameters = sel_inst.pcell_parameters_by_name()
  for k,v in pcell_parameters.items():
    if type(v) in [pya.Path, pya.DPath]:
      return v
      
  return None
  
def set_path_pcell_path(sel_obj, path):
  '''
  Set the <pya.Path> or <pya.DPath> guiding shape in a :sel_obj: <pya.ObjectInstPath> of a PCell
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()

  #assume sel_obj is known to be a pcell
  sel_inst = sel_obj.inst()
  prev_idx = sel_inst.cell.cell_index()
  

  
  pcell_parameters = sel_inst.pcell_parameters_by_name()
  for k,v in pcell_parameters.items():
    if type(v) is type(path):
      sel_inst.change_pcell_parameter(k, path)
      
      #if changes to cell, prune the previous cell      
      if not prev_idx == sel_inst.cell.cell_index(): 
        prev_cell = ly.cell(prev_idx)
        prev_cell.prune_cell()#delete()
        #prev_cell.delete() #this works better?

      return True
          
  return False    
  

def refresh_pcell_path(sel_obj):
  '''
    Refresh a Path guiding shape PCell :sel_obj: <pya.ObjectInstPath> after changes to its guiding shape,
    handling all required clean-up
    '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()

  #update pcell, if valid
  if len(sel_obj.path)>0:
    sel_inst = sel_obj.path[-1].inst()
    if sel_inst.cell.is_pcell_variant():#sel_inst.cell.name in pcell_names:
      prev_idx = sel_inst.cell.cell_index()
      #print('--')
      #print("shape cell idx", prev_cell.cell_index())
      #print("sel inst cell idx", sel_inst.cell.cell_index()) 
      #is a valid pcell, update the appropriate pcell parameter
      
      if sel_inst.cell.name in ['ArrayPath', 'ClearancePath']:
        sel_inst.change_pcell_parameter('s', pya.DPath(sel_obj.shape.path*ly.dbu))
      elif sel_inst.cell.name == 'ROUND_PATH':
        sel_inst.change_pcell_parameter('path', pya.DPath(sel_obj.shape.path*ly.dbu))
      else: #general handling
        pcell_parameters = sel_inst.pcell_parameters_by_name()
        for k,v in pcell_parameters.items():
          if type(v) is pya.Path:
            sel_inst.change_pcell_parameter(k, pya.Path(sel_obj.shape.path))
          elif type(v) is pya.DPath:
            sel_inst.change_pcell_parameter(k, pya.DPath(sel_obj.shape.path*ly.dbu))
        
      #what is the updated cell cell index?
      #print("sel inst cell idx after", sel_inst.cell.cell_index())
      #need to clean up the previous cell
      #prev_cell.delete()
      if not prev_idx == sel_inst.cell.cell_index(): #do not delete cell if no changes made
        prev_cell = ly.cell(prev_idx)
        prev_cell.prune_cell()#delete()



####
# Roundpath Conversions
####

def path2roundpath(sel_obj, verbose=False):
  '''
  Convert a path :sel_obj: <pya.ObjectInstPath> into a ROUND_PATH
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  top_cell = lv.active_cellview().cell

  # check the selected object is a path
  if not sel_obj.shape.is_path() and verbose:
    pya.MessageBox.warning("Selection not a path",
                           "Selection must be a path", pya.MessageBox.Ok)
    return None

  sel_path = sel_obj.shape.path  # this cannot directly be passed into the PCell, need to instantiate a new Path Object
  # ROUND_PATH only accepts DPath (but the coordinates need to be in micron). Path, and Point will not work, the PCell won't show up!
  # https://www.klayout.de/forum/discussion/737/round-path-pcell-in-python-create-cell
  path = pya.DPath()
  path.points = [pya.DPoint(pp.x * ly.dbu, pp.y * ly.dbu) for pp in
                 sel_path.each_point()]  # need to get individual values to avoid rounding
  path.width = sel_path.width * ly.dbu
  layer_info = ly.get_info(sel_obj.layer)

  # remove the selected object
  sel_obj.shape.delete()

  # locate the Basic library
  lib = pya.Library.library_by_name("Basic")
  if lib == None:
    raise Exception("Unknown lib 'Basic'")

  # locate the declaration
  pcell_decl = lib.layout().pcell_declaration("ROUND_PATH");
  if pcell_decl == None:
    raise Exception("Unknown PCell 'ROUND_PATH'")

  # create the PCell variant parameters
  param = {
    "radius": float(pya.Application.instance().get_config("bend-radius")),  # default radius from configs
    "layer": layer_info,  # pya.LayerInfo(2, 0),
    "npoints": int(pya.Application.instance().get_config("bend-points")),
    'path': path
  }

  pv = []
  for p in pcell_decl.get_parameters():
    # print(p.name)
    if p.name in param:
      pv.append(param[p.name])
    else:
      pv.append(p.default)
      # print(p.default)

  # create and place the PCell
  wg_cell_index = ly.add_pcell_variant(lib, pcell_decl.id(), pv)

  t = pya.Trans(pya.Trans.R0, 0, 0)
  inst = top_cell.insert(pya.CellInstArray(wg_cell_index, t))
  ly.cell(wg_cell_index).refresh()
  lv.add_missing_layers()  # add missing layers if an

  # select the new waveguide
  sel_obj = pya.ObjectInstPath()
  sel_obj.top = top_cell.cell_index()
  ie = pya.InstElement(inst)
  sel_obj.append_path(ie)

  # pass the selected object out: selecting here create internal problems with references
  # since we continue to modify other parts of the previous selection
  # lv.select_object(sel_obj)
  return sel_obj


def roundpath2path(sel_obj, verbose=False):
  '''
  Convert a ROUND_PATH :sel_obj: <pya.ObjectInstPath> into a pya.Path
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell

  # get the selected path
  #  object_selection = lv.object_selection
  #  if len(object_selection)<1:
  #    pya.MessageBox.warning("No selection",
  #                                 "Select a path to convert into a round-path", pya.MessageBox.Ok)

  # series of checks ensuring we are selecting an instance of a round path
  # sel_obj = object_selection[0]
  if not sel_obj.is_cell_inst() and verbose:
    pya.MessageBox.warning("Not Round-Path",
                           "Select a round path instance to convert into a path. Hint: Instance frame should be highlighted.",
                           pya.MessageBox.Ok)
    return None, None
  inst = sel_obj.inst()
  if not inst.is_pcell() and verbose:
    pya.MessageBox.warning("Not Round-Path",
                           "Select a round path instant to convert into a path. Hint: Instance frame should be highlighted.",
                           pya.MessageBox.Ok)
    return None, None
  pcell_decl = inst.pcell_declaration()
  if not pcell_decl.name() == "ROUND_PATH" and verbose:
    pya.MessageBox.warning("Not Round-Path",
                           "Select a round path instant to convert into a path. Hint: Instance frame should be highlighted.",
                           pya.MessageBox.Ok)
    return None, None

  # for p in inst.pcell_parameters_by_name():
  #  print(p)

  layer_info = inst.pcell_parameter("layer")
  path = inst.pcell_parameter("path")
  layer_idx = ly.find_layer(layer_info)

  # PCell uses DPoints, but the actual units are Points (very confusing)
  path_points = [pya.DPoint(pp.x, pp.y) for pp in path.each_point()]
  new_path = pya.DPath()
  new_path.points = path_points
  new_path.width = path.width
  shape = cell.shapes(layer_idx).insert(new_path)
  trans = sel_obj.inst().dcplx_trans
  shape.transform(trans)
  # select shape

  del_cell_index = sel_obj.inst().cell.cell_index()
  #  ly.cell(cell_index).delete() #can't delete here, or else we get reference errors
  # sel_obj.inst().delete()

  sel_obj = pya.ObjectInstPath()
  sel_obj.shape = shape
  sel_obj.layer = shape.layer
  sel_obj.top = shape.cell.cell_index()

  # lv.select_object(sel_obj)
  return sel_obj, del_cell_index


def instance_bend(angle):
  '''
  Instance a circular bend PCell with angle :angle:
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  top_cell = lv.active_cellview().cell

  layer_info = ly.get_info(
    ly.find_layer(lv.current_layer.current().source_layer, lv.current_layer.current().source_datatype))

  # locate the Basic library
  lib = pya.Library.library_by_name("Basic")
  if lib == None:
    raise Exception("Unknown lib 'Basic'")

  # locate the declaration
  pcell_decl = lib.layout().pcell_declaration("ARC");
  if pcell_decl == None:
    raise Exception("Unknown PCell 'ARC'")

  radius = float(pya.Application.instance().get_config("bend-radius"))
  wgwidth = float(pya.Application.instance().get_config("wg-width"))
  dbu = ly.dbu

  r1 = radius - wgwidth / 2
  r2 = radius + wgwidth / 2.

  a1 = 0
  a2 = angle  # to be set
  handle1 = pya.Point(r1 * math.cos(a1 * math.pi / 180.) / dbu,
                      r1 * math.sin(a1 * math.pi / 180.) / dbu)  # set in database units
  handle2 = pya.Point(r2 * math.cos(a2 * math.pi / 180.) / dbu, r2 * math.sin(a2 * math.pi / 180.) / dbu)

  # create the PCell variant parameters
  param = {
    "a1": a1,
    "a2": a2,
    "actual_radius1": r1,
    "actual_radius2": r2,
    "radius1": r1,
    "radius2": r2,
    "handle1": handle1,
    "handle2": handle2,
    "actual_start_angle": a1,
    "actual_end_angle": a2,
    "actual_handle1": handle1,
    "actual_handle2": handle2,
    "layer": layer_info,  # pya.LayerInfo(2, 0),
    "npoints": int(pya.Application.instance().get_config("bend-points"))
  }

  pv = []
  for p in pcell_decl.get_parameters():
    print(p.name)
    if p.name in param:
      pv.append(param[p.name])
    else:
      pv.append(p.default)
      # print(p.default)

  # create and place the PCell
  wg_cell_index = ly.add_pcell_variant(lib, pcell_decl.id(), pv)

  # current center
  c = lv.box().center()

  t = pya.Trans(pya.Trans.R0, c.x / dbu, c.y / dbu)  # translate the instance to the current center of view
  inst = top_cell.insert(pya.CellInstArray(wg_cell_index, t))
  ly.cell(wg_cell_index).refresh()
  lv.add_missing_layers()  # add missing layers if any

  # select the new waveguide
  sel_obj = pya.ObjectInstPath()
  sel_obj.top = top_cell.cell_index()
  ie = pya.InstElement(inst)
  sel_obj.append_path(ie)

  # pass the selected object out: selecting here create internal problems with references
  # since we continue to modify other parts of the previous selection
  lv.select_object(sel_obj)
  return sel_obj


def instance_text(text, x=None, y=None):
  '''
  Instance a Text PCell with <string> :text: at position (:x:, :y:)
  '''
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  top_cell = lv.active_cellview().cell

  #  text = "2"
  mag = int(pya.Application.instance().get_config("annotation-mag"))

  # get the current center
  current_box = lv.box()
  c = current_box.center()
  dbu = ly.dbu

  if x is None:
    x = int(round(c.x / dbu))
  if y is None:
    y = int(round(c.y / dbu))

  layer_info = ly.get_info(
    ly.find_layer(lv.current_layer.current().source_layer, lv.current_layer.current().source_datatype))

  # locate the Basic library
  lib = pya.Library.library_by_name("Basic")
  if lib == None:
    raise Exception("Unknown lib 'Basic'")

  # locate the declaration
  pcell_decl = lib.layout().pcell_declaration("TEXT");
  if pcell_decl == None:
    raise Exception("Unknown PCell 'TEXT'")

  # create the PCell variant parameters
  param = {
    "text": text,  # default radius from configs
    "layer": layer_info,  # pya.LayerInfo(2, 0),
    "mag": mag
  }

  pv = []
  for p in pcell_decl.get_parameters():
    #print(p.name)
    if p.name in param:
      pv.append(param[p.name])
      #print(param[p.name])
    else:
      pv.append(p.default)
      #print(p.default)

  # create and place the PCell
  text_cell_index = ly.add_pcell_variant(lib, pcell_decl.id(), pv)

  t = pya.Trans(pya.Trans.R0, x, y)
  inst = top_cell.insert(pya.CellInstArray(text_cell_index, t))
  ly.cell(text_cell_index).refresh()
  lv.add_missing_layers()  # add missing layers if an
