<div id="top"></div>

<br />
<div align="center">
  <h1 align="center">Parallel and Distributed Systems Assignment 2</h1>
  <h3 align="center">Aristotle University of Thessaloniki</h3>
  <h4 align="center">School of Electrical & Computer Engineering</h4>
  <p align="center">
    Contributors: Koro Erika, Vavoulidis Miltiadis-Zacharias
    <br />
    Winter Semester 2021 - 2022
    <br />
    <br />
    <br />
    <br />
  </p>
</div>

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-this-project">About this project</a></li>
    <li><a href="#getting-started">Getting started</a></li>
    <li><a href="#dependencies">Dependencies</a></li>
    <li><a href="#compile-and-run">Complile and run</a></li>
    <li><a href="#data">Data</a></li>
  </ol>
</details>

## About this project

<p align="justify">
  N:Number of points

  p:Number of processes
  
  The objective of this assignment is to re-distribute the points in the processes so as to have in the first p/2 processes, points with distances from the pivot, that are smaller than the median distance and in second half p/2 processes, points with distances bigger than the median distance. We assume that we have number of processes and points that are power of 2. Each process in the end of the program has as many points as it had after the first distribution, N/p.
<br/>
<br/>
  For this exercise we use MPI.
<br/>
<br/>
</p>

## Getting started

To setup this repository on your local machine run the following commands on the terminal:

```console
git clone https://github.com/ErikaKoro/pds-exercise-2.git
```

Or alternatively [*download*](https://github.com/ErikaKoro/pds-exercise-2/archive/refs/heads/main.zip) and extract the zip file of the repository
<br/>
<br/>

## Dependencies
#### 1. Make or Cmake

This project uses make utilities to build and run the executables. 

#### 2. OpenMPI


## Compile and run

### Linux
Simply run `make mpi_sort_build` and after that `make mpi_sort_run`.

To change the number of processes modify the hosts and give to the slots as many processes as you want, considering that the number should be a power of 2.
<br/>

