#
# aroma_harness.py
# 
# Test suite driver
#

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
      print("\t Time to run test: " +  repr(teestTime))
      print("\t Result: \t" + testClass.getResult())
      print("\t Expected: \t" + testClass.getExpectedResult())
      print("\t Test passed: \t" + testClass.testPassed())
      print("="*50)

      self.totalTime += testTime

      if (testClass.testPassed()):
         self.numberOfTestsPassed += 1

   def summary(self):
      print("Total Number of Tests: " + repr(self.numberOfTests))
      print("Total Number of Tests Passed: " + repr(self.numberOfTestsPassed)) 
      print("Total Time: " + repr(self.totalTime))

if __name__=="__main__":
   harness = Harness()


   harness.summary()
  
