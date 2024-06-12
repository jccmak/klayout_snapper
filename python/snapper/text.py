'''
ports.py

snapper v 0.0.2
Copyright 2024 Jason Mak

This file is part of snapper.

snapper is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

snapper is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with snapper. If not, see <https://www.gnu.org/licenses/>.

'''

import pya

def cycle_justification(align, direction_label, cycle):
  lv = pya.Application.instance().main_window().current_view()
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell
   
  object_selection = list(lv.object_selection) #important to ensure ordering stays as expected, or else may behave like a set?
  object_selection = sorted(object_selection, key =lambda x: x.seq)
  
  #cycle = {repr(pya.Text.VAlignTop):pya.Text.VAlignTop,  repr(pya.Text.VAlignCenter):pya.Text.VAlignTop,  repr(pya.Text.VAlignBottom):pya.Text.VAlignCenter,  repr(pya.Text.NoVAlign):pya.Text.VAlignBottom}
  
  lv.transaction("cycle label %s %s"%(align, direction_label))
  
  
  for sel_obj in object_selection:
    if not sel_obj.is_cell_inst(): #is a shape of some sort
      shape = sel_obj.shape
      print(shape)
      if shape.is_text(): #is text
  #        set the justification
        text = shape.text
        cycle_align = cycle[repr(getattr(text, align))] #can't use the index object directly as a key for a dict
        setattr(text, align, cycle_align)
        shape.text = text
  lv.commit()
  
  
def set_magnification(mag):
  
  lv = pya.Application.instance().main_window().current_view() #something weird happens if you don't pass lv from outside. Cell hierarchy selection gets reduced to only what you click directly
  ly = lv.active_cellview().layout()
  cell = lv.active_cellview().cell            
  
  
  
  ##Replace strings of selected cells within the cell hierarchy
  cell_paths = lv.selected_cells_paths(lv.active_cellview().index())
  selected_cell_indices = [cp[-1] for cp in cell_paths]
  
  object_selection = list(lv.object_selection) #important to ensure ordering stays as expected, or else may behave like a set?
  object_selection = sorted(object_selection, key =lambda x: x.seq)
  lv.transaction("set magnification %0.2f"%mag)
  ### If selection with the active layout is not empty
  for sel_obj in object_selection:
    if sel_obj.is_cell_inst():
      inst = sel_obj.inst()
      dcplx_trans = inst.dcplx_trans
      dcplx_trans.mag = mag
      inst.dcplx_trans = dcplx_trans
      #set magnification
    else:
      shape = sel_obj.shape
      #print(shape)
      if shape.is_text():
        text = shape.text
        text.size = int(round(mag/ly.dbu)) #in dbu
        shape.text = text
        
  lv.commit()