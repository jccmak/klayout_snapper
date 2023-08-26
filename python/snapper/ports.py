'''
ports.py

snapper v 0.0.1
Copyright 2023 Jason Mak

This file is part of snapper.

snapper is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

snapper is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with snapper. If not, see <https://www.gnu.org/licenses/>.

'''

import pya
# Enter your Python code here

def gdslayer():
  '''
  Returns the currently configured port gdslayer
  '''
  port_layer = int(round(float(pya.Application.instance().get_config('port-layer'))))
  return port_layer

def get_instance_port_edges(sel_obj):
  '''
  Evaluates the positions of the interface facets marked by the ports on :sel_obj:, a selected instance

  Arguments
  --------
  :sel_obj: <pya.ObjectInstPath> A selected object, assumed to be a cell instance with ports

  Returns
  -------
  :edges: <list> of <pya.Edge> A list of edges marked as a port

  '''
  #points = []
  edges  = []
  if sel_obj.is_cell_inst():
    inst = sel_obj.inst()
    trans = inst.trans #get the transform of the instance
    cell = inst.cell
    port_indices, port_datatypes = ports.occupied_ports(cell)#, sel_obj.layout())
    
    #list out all the edges and points
    for ii in port_indices:
      for ss in cell.shapes(ii).each():
        port_poly = ss.polygon
        port_poly.transform(trans) #transform the poly to the current top level
        
        center, edge = ports.resolve_port(port_poly)
        #print(center, edge)
        #points.append(center)
        edges.append(edge)
  return edges
def resolve_port(port_polygon, direction = False):
  '''
  Given a polygon assumed to be a port, get the edge (and direction, optionally) of the facet marked with the port

  Arguments
  --------
  :port_polygon: <pya.Polygon> the polygon representing the port

  kwargs
  :direction: <bool> Return the direction of the port, additionally

  Returns
  -------
  :point: <pya.Point> the center point of the edge marked as a port
  :edge: <pya.Edge> the marked as a port

  if :direction:

  :direction: <pya.Vector> direction of the input into facet

  '''
  #Ports are defined by an isosceles triangle, with the primary facet being the base (shortest edge)
  point = None
  edge = None
  length = None
  for ee in port_polygon.each_edge():
    if length is None:
      length = ee.length()
      edge = ee
    
    if ee.length()<length:
      length = ee.length()
      edge = ee
  point = (edge.p1 + edge.p2)/2.
  
  if direction is False:
    return point, edge
    
  #calculate the direction as well
  opposing_corner = None
  for pp in port_polygon.each_point_hull():
    if not((pp - edge.p1).length() == 0 or (pp - edge.p1).length() == 0):
      opposing_corner = pp
  direction = opposing_corner - point
  return point, edge, direction

  
  

def ports_info(cell):
  '''
  Given a cell, summarize information about its ports

  Arguments
  --------
  :cell: <pya.Cell> the cell to be evaluated for ports

  Returns
  -------
  :port_indices: <list> of <int> layer index which are port layers
  :port_datatypes: <list> of <int> datatypes for corresponding port layers. Note: all port layers have gdslayer ports.gdslayer()
  :port_empty: <list> of <bool> marks which :port_indices: in the corresponding position are empty
  '''
  ly = cell.layout()
  
  port_layer = int(round(float(pya.Application.instance().get_config('port-layer'))))
  
  port_datatypes = []
  port_indices = [] 
  port_empty = []
  
  layer_indices = ly.layer_indexes() #get all layer index
  for ll in layer_indices:
    info = ly.get_info(ll)
    if info.layer == port_layer:
      port_indices.append(ll)
      port_datatypes.append(info.datatype)
      port_empty.append(cell.shapes(ll).is_empty())
  
  #in sorted order, by datatype
  idx_sorted = sorted(range(len(port_datatypes)), key = lambda x: port_datatypes[x])
  port_datatypes = [port_datatypes[ii] for ii in idx_sorted]
  port_indices = [port_indices[ii] for ii in idx_sorted]
  port_empty = [port_empty[ii] for ii in idx_sorted]
      
  return port_indices, port_datatypes, port_empty
  
def occupied_ports(cell):
  '''
  Given a cell, summarize which ports are occupied

  Arguments
  --------
  :cell: <pya.Cell> the cell to be evaluated for ports

  Returns
  -------
  :port_indices: <list> of <int> layer index which are port layers that are in :cell:
  :port_datatypes: <list> of <int> datatypes for corresponding port layers. Note: all port layers have gdslayer ports.gdslayer()
  '''
  _port_indices, _port_datatypes, _port_empty = ports_info(cell)
  
  port_indices = []
  port_datatypes = []
  for i, is_empty in enumerate(_port_empty):
    if not is_empty:
      port_indices.append(_port_indices[i])
      port_datatypes.append(_port_datatypes[i])
  return port_indices, port_datatypes
  
  


def next_available_port():
  '''
  Return the next unoccupied port in the current cell

  Arguments
  --------

  Returns
  -------
  :port_index_next: <int> layer index of the next unoccupied port :cell: (newly created if not previously available)
  :port_datatype_next: <int> datatype of next unoccupied port
  '''
  lv = pya.Application.instance().main_window().current_view()
  port_layer = int(round(float(pya.Application.instance().get_config('port-layer'))))
  
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell
  
  port_indices, port_datatypes, port_empty = ports_info(cell)
  
  #populate on the next lowest port layer
  if sum(port_empty): #if there are empty layers
    port_empty_min_idx = max(range(len(port_empty)), key=lambda x : port_empty[x]) #argmax
    
    port_index_next = port_indices[port_empty_min_idx]
    port_datatype_next = port_datatypes[port_empty_min_idx]
        
  else: #if all ports are populated, then iterate to next port layer
    if not port_datatypes:
      max_port_datatypes = 0
    else:
      max_port_datatypes = max(port_datatypes)
    
    #create layer for the next datatype
    port_datatype_next = max_port_datatypes + 1
    port_index_next = ly.layer(port_layer,port_datatype_next)
    lv.add_missing_layers()
    
  return port_index_next, port_datatype_next
  
def layer_indices():
  '''
  Returns layer indices which correspond to port layers
  '''
  lv = pya.Application.instance().main_window().current_view()
  
  results = []
  ports_gdslayer = gdslayer()
  for info in lv.each_layer():
    if info.source_layer == ports_gdslayer:
      results.append(info.layer_index())
  return results