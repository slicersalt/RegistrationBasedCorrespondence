// #include "itkTranslationTransform.h"
#include "itkMacro.h"
#include "itkBSplineTransform.h"
#include "itkEuclideanDistancePointMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNormalizedCorrelationPointSetToImageMetric.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkLBFGSBOptimizerv4.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkNormalVariateGenerator.h"
#include "itkImageFileWriter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "itkDisplacementFieldTransform.h"

#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"


#include "vtkDecimatePro.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

#include "MeshToMeshDiffeoRegistrationCLP.h"
 
#include <iostream>
#include <fstream>

namespace 
{
  constexpr unsigned int Dimension = 3;
  using ImageType = itk::Image<float, Dimension>;
  using FixedImageType = ImageType;
  using MovingImageType = ImageType;
    
 
  using MaskPixelType = unsigned char;
 
  using MaskImageType = itk::Image<MaskPixelType, Dimension>;
 
  using FixedPointSetType = itk::PointSet<float, Dimension>;
 
  using PixelType = signed short;
 
  
 
  using MaskPixelType = unsigned char;
 
  using MaskImageType = itk::Image<MaskPixelType, Dimension>;

  using TransformType = itk::BSplineTransform<double,3,3>;
 
  using ParametersType = TransformType::ParametersType;
 
 
  using OptimizerType = 
    itk::LBFGSBOptimizerv4;
   
  using LinearInterpolatorType =
    itk::LinearInterpolateImageFunction<ImageType, double>;
 
  using MetricType =
    itk::MeanSquaresImageToImageMetricv4<ImageType,ImageType>;

  using OptimizerScalesType = OptimizerType::ScalesType;
 
  using RegistrationType =
    itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType>;

    using MeshType = itk::Mesh<float, Dimension>;

class CommandIterationUpdate : public itk::Command
{
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<CommandIterationUpdate>;
  itkNewMacro(CommandIterationUpdate);
 
protected:
  CommandIterationUpdate() = default;
 
  using InternalImageType = itk::Image<float, 2>;
  using VectorPixelType = itk::Vector<float, 2>;
  using DisplacementFieldType = itk::Image<VectorPixelType, 2>;
 
  using RegistrationFilterType =
    itk::DiffeomorphicDemonsRegistrationFilter<InternalImageType,
                                                 InternalImageType,
                                                 DisplacementFieldType>;
 
public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }
 
  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    const auto * filter = static_cast<const RegistrationFilterType *>(object);
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    std::cout << filter->GetMetric() << std::endl;
  }
};

ImageType::Pointer meshToImage(MeshType::Pointer mesh)
{
  auto bounds = mesh->GetBoundingBox()->GetBounds();

  ImageType::SpacingType spacing;
  spacing[0] = (bounds[1] - bounds[0]) / 90;
  spacing[1] = (bounds[3] - bounds[2]) / 90;
  spacing[2] = (bounds[5] - bounds[4]) / 90;

  ImageType::PointType origin;
  origin[0] = bounds[0] - 5*spacing[0];
  origin[1] = bounds[2] - 5*spacing[1];
  origin[2] = bounds[4] - 5*spacing[2];
  
  ImageType::SizeType size;
  size[0] = 100;
  size[1] = 100;
  size[2] = 100;

  using MeshToImageType = itk::TriangleMeshToBinaryImageFilter<MeshType,ImageType>;
  auto meshToImage = MeshToImageType::New();
  meshToImage->SetInput(mesh);
  meshToImage->SetOrigin(origin);
  meshToImage->SetSpacing(spacing);
  meshToImage->SetSize(size);
  meshToImage->Update();

  using DistanceType = itk::SignedMaurerDistanceMapImageFilter<ImageType,ImageType>;
  auto distance = DistanceType::New();
  distance->SetInput(meshToImage->GetOutput());
  distance->Update();

  auto image = distance->GetOutput();
  return image;
}

ImageType::Pointer meshToImage(MeshType::Pointer mesh, ImageType::Pointer referenceImage)
{
  using MeshToImageType = itk::TriangleMeshToBinaryImageFilter<MeshType,ImageType>;
  auto meshToImage = MeshToImageType::New();
  meshToImage->SetInput(mesh);
  meshToImage->SetOrigin(referenceImage->GetOrigin());
  meshToImage->SetSpacing(referenceImage->GetSpacing());
  meshToImage->SetSize(referenceImage->GetLargestPossibleRegion().GetSize());
  meshToImage->Update();

  using DistanceType = itk::SignedMaurerDistanceMapImageFilter<ImageType,ImageType>;
  auto distance = DistanceType::New();
  distance->SetInput(meshToImage->GetOutput());
  distance->Update();

  auto image = distance->GetOutput();
  return image;
}

MeshType::Pointer polyDataToMesh(vtkSmartPointer<vtkPolyData> pd)
{
  auto mesh = MeshType::New();

  int numberOfPoints = pd->GetNumberOfPoints();
  mesh->GetPoints()->Reserve( numberOfPoints );

  for (int p = 0; p < numberOfPoints; p++)
  {
    double* point = pd->GetPoint(p);
    mesh->SetPoint( p, MeshType::PointType( point ) );
  }

  int numberOfTriangles = pd->GetNumberOfCells();
  mesh->GetCells()->Reserve( numberOfTriangles );

  typedef MeshType::CellType CellType;
  typedef itk::TriangleCell<CellType> TriangleCellType;

  for (int c = 0; c < numberOfTriangles; c++)
  {
    unsigned long pointIds[3];
    pointIds[0] = pd->GetCell(c)->GetPointIds()->GetId(0);
    pointIds[1] = pd->GetCell(c)->GetPointIds()->GetId(1);
    pointIds[2] = pd->GetCell(c)->GetPointIds()->GetId(2);
    MeshType::CellAutoPointer itkCell;
    TriangleCellType *tcell = new TriangleCellType;
    TriangleCellType::PointIdentifier itkPts[3];
    for (int ii = 0; ii < 3; ++ii)
    {
      itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(pointIds[ii]);
    }
    tcell->SetPointIds(itkPts);
    itkCell.TakeOwnership(tcell);
    mesh->SetCell(c,itkCell);
  }  

  return mesh;
}

int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;
 
  // Read in meshes

  auto templateReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  templateReader->SetFileName(templateMeshFile.c_str());
  templateReader->Update();
  auto templateMesh = templateReader->GetOutput();

  auto targetReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  targetReader->SetFileName(targetMeshFile.c_str());
  targetReader->Update();
  auto targetMesh = targetReader->GetOutput();

  // Convert target to ITK mesh
  std::cout << "Converting to ITK mesh" << std::endl;
  auto targetITKMesh = polyDataToMesh(targetMesh);
  auto templateITKMesh = polyDataToMesh(templateMesh);

  std::cout << "done" << std::endl;

  // Convert meshes to images
  // This is opposite from what you would expect but due 
  // to the way ITK handles transforms this will correctly
  // yield a registration of the template into the target
  // rather than the target into the template
  
  auto fixedImage = meshToImage(templateITKMesh);
  auto movingImage = meshToImage(targetITKMesh,fixedImage);
  
using VectorPixelType = itk::Vector<double, Dimension>;
using DisplacementFieldType = itk::Image<VectorPixelType, Dimension>;
using RegistrationFilterType =
  itk::DiffeomorphicDemonsRegistrationFilter<ImageType,
                                                ImageType,
                                                DisplacementFieldType>;
  RegistrationFilterType::Pointer filter = RegistrationFilterType::New();

 
 CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  filter->AddObserver(itk::IterationEvent(), observer);
 
  filter->SetFixedImage(fixedImage);
  filter->SetMovingImage(movingImage);
 
  filter->SetNumberOfIterations(200);
  filter->SetStandardDeviations(1.0);
  filter->Update();
  
  std::cout << "done" << std::endl;

  using DisplacementFieldTransformType =
    itk::DisplacementFieldTransform<double, Dimension>;
  auto transform = DisplacementFieldTransformType::New();
  transform->SetDisplacementField(filter->GetOutput());

  using PointIdentifierType = MeshType::PointIdentifier;
  MeshType::PointType transformedPoint;
  for( PointIdentifierType pointId = 0; pointId < templateMesh->GetNumberOfPoints(); ++pointId )
  {
    transformedPoint[0] = templateMesh->GetPoint(pointId)[0];
    transformedPoint[1] = templateMesh->GetPoint(pointId)[1];
    transformedPoint[2] = templateMesh->GetPoint(pointId)[2];
    transformedPoint = transform->TransformPoint( transformedPoint );
    templateMesh->GetPoints()->SetPoint(pointId, transformedPoint[0], transformedPoint[1], transformedPoint[2]);
  }

  
  auto meshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  meshWriter->SetInputData(templateMesh);
  meshWriter->SetFileName(registeredTemplateFile.c_str());
  meshWriter->Update();

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  try
    {
    // This filter handles all types on input, but only produces
    // signed types
    DoIt(argc, argv);
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
