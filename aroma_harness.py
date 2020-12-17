#
# aroma_harness.py
# 
# Test suite driver
#

import time
import fpformat

import aroma
import aroma_molecule

class TestRun:
  def runAndGetTime(self):
      pass

  def getName(self):
      pass

  def getResult(self):
      pass

  def getExpectedResult(self):
      pass

  def testPassed(self):
      pass

class SimpleAromaRunZScan:
   def __init__(self):
      self.result = []
      self.expectedResult = []

   def getName(self):
      return "SimpleAromRun - z-scan"

   def runAndGetTime(self):
      start = time.time()

      aroma.aroma("input/benz")

      self.result = open("output/benz.armlog").readlines()
      self.expectedResult = open("testout/benz.armlog").readlines()

      end = time.time()

      return (end-start)

   def getExpectedResult(self):
      return repr(self.expectedResult)

   def getResult(self):
      return repr(self.result)

   def testPassed(self):
      return self.result == self.expectedResult

class ZMatToCartesian:
   def __init__(self):
      self.result = []
      self.expectedResult = []

   def getName(self):
      return "ZMatToCartesian"

   def runAndGetTime(self):
      start = time.time()

      inpfl = open("input/zmat.in")
      zmatLines = list(map(lambda x: x.strip().split(), inpfl.readlines()[5:]))
      inpfl.close()
      
      print(zmatLines)

      aroma_molecule.generateCartesianFromZmat(zmatLines)
   

      end = time.time()

      return (end-start)

   def getExpectedResult(self):
      return repr(self.expectedResult)

   def getResult(self):
      return repr(self.result)

   def testPassed(self):
      return self.result == self.expectedResult
      
class Harness:
   def __init__(self):
      self.numberOfTests = 0
      self.numberOfTestsPassed = 0
      self.totalTime = 0.0

   def run(self, testClass):
      self.numberOfTests += 1

      print("="*50)
      print("Running Test: " + testClass.getName())
      testTime = testClass.runAndGetTime()
      print("\t Time to run test: " +  fpformat.fix(testTime, 4))
      print("\t Result: \t" + testClass.getResult())
      print("\t Expected: \t" + testClass.getExpectedResult())
      print("\t Test passed: \t" + repr(testClass.testPassed()))
      print("="*50)

      self.totalTime += testTime

      if (testClass.testPassed()):
         self.numberOfTestsPassed += 1

   def summary(self):
      print("Total Number of Tests: " + repr(self.numberOfTests))
      print("Total Number of Tests Passed: " + repr(self.numberOfTestsPassed)) 
      print("Total Time: " + fpformat.fix(self.totalTime, 4))

if __name__=="__main__":
   harness = Harness()

   harness.run(SimpleAromaRunZScan())
   harness.run(ZMatToCartesian())

   harness.summary()
  
