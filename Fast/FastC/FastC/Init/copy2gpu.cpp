/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
//# include <numa.h>
//# include </usr/src/kernels/2.6.32-573.12.1.el6.x86_64/include/linux/numa.h>

# include "FastC/fastc.h"
//# include "converter.h"
# include "kcore.h"
# include "FastC/param_solver.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Init data in parallel openmp to improve data placement on numa machine */
//=============================================================================
PyObject* K_FASTC::copy2gpu(PyObject* self, PyObject* args)
{
  PyObject* targetArray;  E_Int nvar;  E_Int size; 

#ifdef E_DOUBLEINT 
  if (!PyArg_ParseTuple(args, "Oll", &targetArray, &size, &nvar ))
#else
  if (!PyArg_ParseTuple(args, "Oii", &targetArray, &size, &nvar ))
#endif
  {
    return NULL;
  }

  FldArrayF* target;
  K_NUMPY::getFromNumpyArray(targetArray, target); E_Float* ipttarget = target->begin();

#ifdef E_GPU
E_Int sizetab= nvar*size;
//#pragma omp target map(to ipttarget[:sizetab])
printf("Copy2GPU dans compact\n");
#endif

  RELEASESHAREDN( targetArray  , target  );

  Py_INCREF(Py_None);
  return Py_None;

}
