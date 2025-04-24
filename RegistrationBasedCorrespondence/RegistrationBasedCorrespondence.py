from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import vtk
import qt
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule, ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleWidget, ScriptedLoadableModuleTest,
)


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
    self.parent.contributors = ["David Allemang (Kitware), Beatriz Paniagua (Kitware)"]
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
    self.logic: Optional[RegistrationBasedCorrespondenceLogic] = None

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer)
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/RegistrationBasedCorrespondence.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    self.methodButtonGroup = qt.QButtonGroup()
    self.methodButtonGroup.setExclusive(True)
    self.methodButtonGroup.addButton(self.ui.DiffeoButton)
    self.methodButtonGroup.addButton(self.ui.BSplineButton)

    self.logic = RegistrationBasedCorrespondenceLogic()

    self.ui.ApplyButton.connect('clicked(bool)', self.onApplyButton)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    try:
      self.logic.run(
        Path(self.ui.TemplateMesh.currentPath),
        Path(self.ui.InputDirectory.directory),
        Path(self.ui.OutputDirectory.directory),
        self.methodButtonGroup,
        self.ui.IterationsSlider.value,
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

  def run(self, template: Path, data: Path, output: Path, methodButtonGroup, iterations):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param template:
    :param data:
    :param output:
    """
    if not template.exists():
      raise ValueError('Template file is not valid.')

    if not data.is_dir():
      raise ValueError('data directory is not valid.')

    if output.exists() and not output.is_dir():
      raise ValueError('output directory is not valid.')

    # Choose registration method
    method = methodButtonGroup.checkedButton().objectName

    if method == 'DiffeoButton':
      cliToRun = slicer.modules.meshtomeshdiffeoregistration
    else:
      cliToRun = slicer.modules.meshtomeshregistration

    logging.info('Processing started')

    storage_check = slicer.vtkMRMLModelStorageNode()

    files = os.listdir(data)
    results = [template]
    for file in files:
      if not storage_check.GetSupportedFileExtension(file):
        logging.info("Skipping file %r. Invalid model format.", file)
        continue

      inputFile = os.path.join(data,file)
      outputFile = os.path.join(output,file)

      cliParams = {
          'templateMeshFile': os.fspath(template),
          'targetMeshFile': inputFile,
          'registeredTemplateFile': outputFile,
          'iterations': iterations,
      }
      cliNode = slicer.cli.runSync(cliToRun, None, cliParams)

      results.append(outputFile)

    # Load shapes into SPV
    slicer.modules.shapepopulationviewer.widgetRepresentation().deleteModels()
    for result in results:
      _, ext = os.path.splitext(result)
      pdr = vtk.vtkXMLPolyDataReader() if ext == ".vtp" else vtk.vtkPolyDataReader()
      pdr.SetFileName(result)
      pdr.Update()

      # Add scalar fields for QC in SPV
      scalars = vtk.vtkIntArray()
      scalars.SetName("QC_RegistrationBasedCorrespondence")  
      for i in range(pdr.GetOutput().GetNumberOfPoints()):
          scalars.InsertNextValue(i)  

      pdr.GetOutput().GetPointData().SetScalars(scalars)
      pdr.Update()

      slicer.modules.shapepopulationviewer.widgetRepresentation().loadModel(pdr.GetOutput(),result)

    slicer.util.selectModule(slicer.modules.shapepopulationviewer)

    logging.info('Processing completed')


#
# RegistrationBasedCorrespondenceTest
#

class RegistrationBasedCorrespondenceTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_RegistrationBasedCorrespondence1()

  def test_RegistrationBasedCorrespondence1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    import tempfile
    import os

    logic = RegistrationBasedCorrespondenceLogic()

    with tempfile.TemporaryDirectory() as tempdir:
      tempdir = Path(tempdir)

      content1 = os.urandom(32)
      content2 = os.urandom(32)

      data = tempdir / 'data'
      data.mkdir()
      (data / 'file').write_bytes(content1)
      (data / 'sub').mkdir()
      (data / 'sub' / 'file').write_bytes(content2)

      output = tempdir / 'output'

      logic.run(data, output)

      self.assertTrue(output.exists())
      self.assertTrue((output / 'file').exists())
      self.assertEqual((output / 'file').read_bytes(), content1)

      self.assertTrue((output / 'sub').exists())
      self.assertTrue((output / 'file').exists())
      self.assertEqual((output / 'sub' / 'file').read_bytes(), content2)

    self.delayDisplay('Test passed')
