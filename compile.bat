rem /*
rem  * Copyright (c) 2011-2012 Vikram Sundar.
rem  * All Rights Reserved.
rem  */
rmdir /s /q org
del /f ProteinDock.jar
javac -d . src\java\org\vikramdock\ProteinDockPredict.java src\java\org\vikramdock\Atom.java src\java\org\vikramdock\ProteinStruct.java src\java\org\vikramdock\TestCase.java src\java\org\vikramdock\Constants.java src\java\org\vikramdock\Bond.java src\java\org\vikramdock\TestCaseGenerator.java src\java\org\vikramdock\TestCaseComparator.java src\java\org\vikramdock\TestCaseStore.java
jar cvf ProteinDock.jar org
