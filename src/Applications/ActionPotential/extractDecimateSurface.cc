#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;

int main(int argc, char* argv[])
{
  if( argc != 3 ) {
      cout << "Usage: indent modelName decimationTarget" << endl;
      return(0);
  }

  string modelName = argv[1];

  string inputFileName = modelName + ".vtk";

  // read in file
  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );

  // extract boundary surface
  vtkGeometryFilter * geometry = vtkGeometryFilter::New();
  geometry->SetInput( reader->GetOutput() );

  // send through normals filter to ensure that triangle orientations
  // are consistent
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  normals->SetInput( geometry->GetOutput() );
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();

  // get the mesh
  vtkPolyData * mesh = normals->GetOutput();
  mesh->Update();
  
  std::cout << "Surface Mesh: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // convert to triangle mesh if necessary
  vtkTriangleFilter * triangles = vtkTriangleFilter::New();
  triangles->SetInput( normals->GetOutput() );
  triangles->Update();

  mesh = triangles->GetOutput();
  mesh->Update();
  std::cout << "Surface Triangulation: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // write basic triangle mesh to a file
  vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
  writer->SetInput( triangles->GetOutput() );
  
  string outputFileName = modelName + "-surface.vtk";
  writer->SetFileName( outputFileName.c_str() );
//   writer->SetFileTypeToBinary();
  writer->SetFileTypeToASCII();
  writer->Write();

  // decimate using vtkDecimatePro
  //
  //  http://www.vtk.org/doc/release/5.4/html/a00397.html
  //
  // Note: two other similar decimation filters are 
  //  vtkQuadricClustering: http://www.vtk.org/doc/release/5.4/html/a01444.html
  //  vtkQuadricDecimation: http://www.vtk.org/doc/release/5.4/html/a01446.html
  //
  vtkDecimatePro * decimate = vtkDecimatePro::New();
  decimate->SetInput( triangles->GetOutput() );
  // How much to decimate, i.e., 0.9 tries to remove 90% of the triangles
  double target = atof(argv[2]);
  decimate->SetTargetReduction( target );
  decimate->PreserveTopologyOn();
  // check vtk online documentation for many many other options for
  // this filter
  decimate->Update();
  
  mesh = decimate->GetOutput();
  mesh->Update();
  std::cout << target*100.0 << "\% Decimated Triangulation: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // write decimated mesh to file
  writer->SetInput( decimate->GetOutput() );
  char tmp[20];
  sprintf(tmp,"%.2g",target);
  outputFileName = modelName + "-decimated-" + tmp + ".vtk";
  writer->SetFileName( outputFileName.c_str() );
  writer->Write();

  return 0;
}

