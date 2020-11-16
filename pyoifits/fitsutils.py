from astropy.io import fits
import numpy as np
import re 
from collections import OrderedDict 

column_dtype = {
    'B': 'B', 'h': 'I', 'i': 'J', 'l': 'K', 
    'f': 'E', 'd': 'D', 'F': 'C', 'D': 'M',
    '?': 'L', 'b': 'L', 'U': 'A', 'S': 'A'
}

def ascolumn(data=None, *, name=None, format=None, unit=None, dim=None,
             null=None, bzero=None, bscale=None, **kwarg):
    """
Return data as a FITS column.

Keyword arguments
-----------------
data:
        Data. Either an array or a FITS column.  If data is None, the format
        must be given.

name (str):
        Name of the column

unit (str):
        Unit

dim (tuple):
        Dimension of column cells.

format (str):
        Column format. If absent or incomplete (letter without dimension)
        it will be completed using data shape and type.

null (int):
        Null value for integer types.

    """

    if unit in ['', 'any']:
        unit = None
     
    # If a FITS column is given, look if we need to transform it

    if isinstance(data, fits.Column):
        
        # Look if it matches specifications
        if ((name is None or data.name == name) and
            (unit is None or data.unit == unit) and
            (null is None or data.null == null) and
            (bzero is None or data.bzero == bscale) and
            (bscale is None or data.bscale == bscale) and
            (dim is None or data.dim == dim)):

            # For format a few cases are possible
            # 1. no format is specified
            # 2. format contains only one letter, no data size
            # 3. format contains data size
            if (format is None or
                len(format) == 1 and data.format[-1] == format or
                data.format == format):
                return data

        # we revert to array and will try to reshape column
        array = data.array
        if name is None:
            name = data.name
        if unit is None:
            unit = data.unit
        if null is None:
            null = data.null
        if dim is None:
            dim = data.dim
    
    else:
        array = data
   
    # If array is None
    if array is None:
        col = fits.Column(name=name, unit=unit, format=format, null=null,
                dim=dim)
        return col

    # Shape is determined from the array 
    array = np.asarray(array)
    shape = array.shape
    size = int(np.prod(shape[1:]))
   
    # Determine TFORM keyword 
    if format is None:
        dtype = array.dtype.char
        format = 'A' if dtype in 'SU' else column_dtype(dtype) 
    
    if len(format) == 1:
        len_ = size
        if format[-1] == 'A':
            array = np.asarray(array, dtype=str)
            len_ *= array.itemsize // array.dtype.alignment
        format = f"{len_}{format}"
    else:
        len_ = size
        if format[-1] == 'A':
            len_ *= int(format[:-1]) // size
        if len_ != int(format[:-1]):
            print(array)
            shape = array.shape[1:]
            txt = f"Shape {shape} does not match format {format} for {name}"
            raise RuntimeError(txt)

    # Determine TDIM keyword. 
    tdim = None
    if format[-1] == 'A' and len(shape) > 1:
        tdim = (len_, *shape[1:])
    elif len(shape) > 2:
        tdim = shape[1:]

    col = fits.Column(array=array, format=format, name=name, 
            unit=unit, dim=dim, null=null)

    return col

def merge_columns(*hdus):

    # Find all column names keeping order of appearance in hdus
    colnames = [col.name for hdu in hdus for col in hdu.columns]
    colnames = list(OrderedDict.fromkeys(colnames))

    # for each column, find all HDUs it is included and keep the
    # columns with largest data size.  
    columns = []
    for name in colnames:
        col = None
        for hdu in hdus:
            if name in hdu.columns.names:
                new_col = hdu.columns[name]
                if col is None:
                    col = new_col
                elif col.array.itemsize < new_col.array.itemsize:
                    col = new_col
        columns.append(col)

    return colnames, columns

def merge_rows(*rows, id_name=None, equality=lambda x,y: x==y):

    if id_name is None:
        return rows, [{}] * len(rows)

    # New unused IDs (in a sequence, avoiding the ones in either
    rows1 = rows[0]
    id1 = rows1[id_name].tolist()
    kept = [rows1]
    maps = [{}]

    for rows2 in rows[1:]:

        id2 = rows2[id_name]
        candidate_id2 = set(range(1, 1 + len(id1) + len(id2)))
        candidate_id2 -= set([*id1, *id2])
        candidate_id2 = sorted(list(candidate_id2))

        index_map = {}
        kept_lines = []
        for j, row2 in enumerate(rows2):

            equal = [equality(row2, row1) for rows1 in kept for row1 in rows1]
            index_equal = np.argwhere(equal)
            if len(index_equal):
                index_map[id2[j]] = id1[index_equal[0,0]]
                continue

            kept_lines.append(j)
            value = candidate_id2.pop(0) if id2[j] in id1 else id2[j]
            index_map[id2[j]] = value
            id1.append(value)

        kept2 = rows2[kept_lines]
        index_map = {o: n for o, n in index_map.items() if o != n}

        kept.append(kept2)
        maps.append(index_map)

    return kept, maps

def merge_fits_headers(*headers, req_keys=[]):

    header = headers[0].copy()

    for h in headers[1:]:
        for card in h.cards:
            (name, value, comment) = card
            if name in ['', 'HISTORY', 'COMMENT']:
                header.append(card)
            elif name not in header:
                if 'HISTORY' in header:
                    header.insert('HISTORY', card)
                else:
                    header.append(card)
            elif value != header[name]:
                if type(header[name]) == str:
                    header[name] = 'MULTI'
                elif name not in req_keys:
                    header.remove(name)

    return header

