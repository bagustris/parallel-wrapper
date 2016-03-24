#!/bin/bash

myid=${OMPI_COMM_WORLD_RANK}

touch file$myid
