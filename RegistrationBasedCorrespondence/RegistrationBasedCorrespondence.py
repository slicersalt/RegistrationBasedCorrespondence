import logging
import shutil
from pathlib import Path

import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


#
# RegistrationBasedCorrespondence
#

class RegistrationBasedCorrespondence(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Registration-based Correspondence"
    self.parent.categories = ["Shape Creation"]
    self.parent.dependencies = []
    self.parent.contributors = ["David Allemang (Kitware)"]
    self.parent.helpText = (
      "Generate Correspondence of consistent non-spherical topologies."
    )
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = (
      "This file was originally developed by Jean-Christophe Fillion-Robin, "
      "Kitware Inc., Andras Lasso, PerkLab, and Steve Pieper, Isomics, Inc. "
      "and was partially funded by NIH grant 3P41RR013218-12S1."
    )


#
# RegistrationBasedCorrespondenceWidget
#

class RegistrationBasedCorrespondenceWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    self.logic = None

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer)
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/RegistrationBasedCorrespondence.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    self.logic = RegistrationBasedCorrespondenceLogic()

    self.ui.ApplyButton.connect('clicked(bool)', self.onApplyButton)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    try:
      self.logic.run(
        Path(self.ui.InputDirectory.directory),
        Path(self.ui.OutputDirectory.directory)
      )
    except Exception as e:
      slicer.util.errorDisplay("Failed to compute results: {}".format(e))
      import traceback
      traceback.print_exc()


#
# RegistrationBasedCorrespondenceLogic
#

class RegistrationBasedCorrespondenceLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self, data: Path, output: Path):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param data:
    :param output:
    """

    if not data.is_dir():
      raise ValueError('data directory is not valid.')

    if output.exists() and not output.is_dir():
      raise ValueError('output directory is not valid.')

    logging.info('Processing started')

    if output.exists():
      shutil.rmtree(output)
    shutil.copytree(str(data), str(output))

    logging.info('Processing completed')

# #
# # RegistrationBasedCorrespondenceTest
# #
#
# class RegistrationBasedCorrespondenceTest(ScriptedLoadableModuleTest):
#   """
#   This is the test case for your scripted module.
#   Uses ScriptedLoadableModuleTest base class, available at:
#   https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
#   """
#
#   def setUp(self):
#     """ Do whatever is needed to reset the state - typically a scene clear will be enough.
#     """
#     slicer.mrmlScene.Clear(0)
#
#   def runTest(self):
#     """Run as few or as many tests as needed here.
#     """
#     self.setUp()
#     self.test_RegistrationBasedCorrespondence1()
#
#   def test_RegistrationBasedCorrespondence1(self):
#     """ Ideally you should have several levels of tests.  At the lowest level
#     tests should exercise the functionality of the logic with different inputs
#     (both valid and invalid).  At higher levels your tests should emulate the
#     way the user would interact with your code and confirm that it still works
#     the way you intended.
#     One of the most important features of the tests is that it should alert other
#     developers when their changes will have an impact on the behavior of your
#     module.  For example, if a developer removes a feature that you depend on,
#     your test should break so they know that the feature is needed.
#     """
#
#     self.delayDisplay("Starting the test")
#
#     # Get/create input data
#
#     import SampleData
#     inputVolume = SampleData.downloadFromURL(
#       nodeNames='MRHead',
#       fileNames='MR-Head.nrrd',
#       uris='https://github.com/Slicer/SlicerTestingData/releases/download/MD5/39b01631b7b38232a220007230624c8e',
#       checksums='MD5:39b01631b7b38232a220007230624c8e')[0]
#     self.delayDisplay('Finished with download and loading')
#
#     inputScalarRange = inputVolume.GetImageData().GetScalarRange()
#     self.assertEqual(inputScalarRange[0], 0)
#     self.assertEqual(inputScalarRange[1], 279)
#
#     outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
#     threshold = 50
#
#     # Test the module logic
#
#     logic = RegistrationBasedCorrespondenceLogic()
#
#     # Test algorithm with non-inverted threshold
#     logic.run(inputVolume, outputVolume, threshold, True)
#     outputScalarRange = outputVolume.GetImageData().GetScalarRange()
#     self.assertEqual(outputScalarRange[0], inputScalarRange[0])
#     self.assertEqual(outputScalarRange[1], threshold)
#
#     # Test algorithm with inverted threshold
#     logic.run(inputVolume, outputVolume, threshold, False)
#     outputScalarRange = outputVolume.GetImageData().GetScalarRange()
#     self.assertEqual(outputScalarRange[0], inputScalarRange[0])
#     self.assertEqual(outputScalarRange[1], inputScalarRange[1])
#
#     self.delayDisplay('Test passed')
