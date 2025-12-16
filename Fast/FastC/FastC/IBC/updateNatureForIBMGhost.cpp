# include "fastc.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* update cellN for IBM (zone, in place) */
// ============================================================================
PyObject* K_FASTC::_updateNatureForIBMGhost(PyObject* self, PyObject* args)
{
  PyObject* zone;
  char* GridCoordinates; char* FlowSolutionNodes;
  char* FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, O_ SSS_,
                    &zone, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
      return NULL;
  }
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(
    zone, 0, 1, varString, fields, locs, im, jm, km, 
    cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
    FlowSolutionNodes, FlowSolutionCenters);

  if ( res != 1 ) 
  {
    if ( res == 2 )
    { 
      delete [] eltType; delete [] varString;
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "updateNatureForIBM: not valid for unstructured grids.");
      return NULL;
    }
    else 
    {
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "updateNatureForIBM: invalid zone.");
      return NULL;
    }
  }  
  E_Int poscellni = K_ARRAY::isNamePresent("cellN", varString);      //cellNIBC au moment de l'appel
  E_Int poscellnc = K_ARRAY::isNamePresent("cellNSave1", varString); //cellNChim      
  if ( poscellni == -1 || poscellnc == -1)
  {
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "updateNatureForIBMGhost: zone must contain cellNIBC and cellNChim variables.");
    return NULL;
  }
  
  E_Float* ptrCellNIBC  = fields[poscellni];
  E_Float* ptrCellNChim = fields[poscellnc];
  E_Int imc = K_FUNC::E_max(1,im-1);
  E_Int jmc = K_FUNC::E_max(1,jm-1);
  E_Int kmc = K_FUNC::E_max(1,km-1);
  E_Int ncells =  imc*jmc*kmc;
#pragma omp parallel default(shared)
  {
#pragma omp for 
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      E_Float& cellNChim = ptrCellNChim[ind];
      E_Float& cellNIBC = ptrCellNIBC[ind];

      E_Float tol     = 1.e-8;

      if ( K_FUNC::fEqualZero(cellNChim-1.,tol) == true )
      {
        if      ( K_FUNC::fEqualZero(cellNIBC    ,tol) == true ) {cellNChim = -3.;}//~ blanked
        else if ( K_FUNC::fEqualZero(cellNIBC -2.,tol) == true ) {cellNChim =  3.;}//not a donor
      }
      else if ( K_FUNC::fEqualZero(cellNChim-2.,tol) == true )
      {
        if     ( K_FUNC::fEqualZero(cellNIBC -2.,tol) == true ) {cellNChim = -3.;}//blanked
        else if( K_FUNC::fEqualZero(cellNIBC    ,tol) == true ) {cellNChim = -3.;}//blanked
      }
      else if ( K_FUNC::fEqualZero(cellNChim ,tol) == true )
      {
        if      ( K_FUNC::fEqualZero(cellNIBC -1.,tol) == true ) {cellNIBC = -3.;}//blanked
      }
      
      // c est commente : on suppose que les corps IBM n intersectent pas les corps Chimere
      // if ( cellNFront != 0.)
      // {
      //   if (cellNIBC == -3.) cellNFront = 0.;
      // }
    }
  }
  delete [] eltType; delete [] varString;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  Py_INCREF(Py_None); 
  return Py_None;
}
