# To run the gp-unit tests, execute them from the gp-unit dir but with
# this gpunit.properties file

cd ../../../util/gpunit/
ant -Dgpunit.properties=../../MutSigCV/gpunit/v1/gpunit.properties gpunit
