rem /*
rem  * Copyright (c) 2011-2012 Vikram Sundar.
rem  * All Rights Reserved.
rem  */
rmdir /s /q org
del /f ProteinDock.jar
javac -d . src\java\org\vikramdock\*.java
jar cvf ProteinDock.jar org
