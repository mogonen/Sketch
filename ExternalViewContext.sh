#!/bin/sh
if [ -z "$HUESPACE3_HOME" ]; then
	HUESPACE_SDK=../../../..
	echo "HUESPACE3_HOME is not set, using $HUESPACE_SDK"
	HUESPACE3_HOME=$HUESPACE_SDK LD_LIBRARY_PATH=$HUESPACE_SDK/lib/ ./ExternalViewContext
else
	echo "HUESPACE3_HOME is set to $HUESPACE3_HOME"
	LD_LIBRARY_PATH=$HUESPACE3_HOME/lib/ ./ExternalViewContext
fi
