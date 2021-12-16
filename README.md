<div id="top"></div>
<br />
<div align="center">
  <h1 align="center">Parallel and Distributed Systems Assignment 2</h1>
  <h3 align="center">Aristotle University of Thessaloniki</h3>
  <h4 align="center">School of Electrical & Computer Engineering</h4>
  <p align="center">
    Contributors: Koro Erika, Vavoulidis Miltiadis - Zacharias
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
    <li><a href="#project-directory-layout">Project directory layout</a></li>
    <li><a href="#compile-and-run">Complile and run</a></li>
    <li><a href="#graphs">Graphs</a></li>
    <li><a href="#useful-links">Useful links</a></li>
  </ol>
</details>
##About this project
<p align="justify">
The objective of this project is to distribute an amount of given points (given their coordinates) by their median to a number of processes using Message Passing Interface (MPI) in C programming language.At first, we were given a number of points (consider this as "N") distributed in some processes (consider number of processes as "P") where every process holds locally N/P points. You may consider the dimensions of the points as "d".The distribution took place using the following algorithm:
<br/>
At first, we chose a process to be the _**master**_ process which, sort of, controls the data/info that we want to be exhanged between all the processes. Then, the master chooses one of its points to be the **pivot** and announces it to all the other processes. Every process calculates the Euclidean distance between each one of their points and the announced pivot and after that the master gathers all this data. The master now, finds the median of all the gathered distances, using the **quickselect** routine. Aftewards, all the median is announced to all the processes and from now on it is up to every process to seperate their respective points to those whose distances are smaller and to those whose distances are greater than the median.The objective of all that stuff is to create an array of exchanges (of points) so that the points will be re-distributed
<br/>
<br/>
After the distribution is done:
<br/>
Note 1: All the processes must locally hold exactly N/P points.
<br/> 
Note 2: The smaller the rank of a process the smaller the distance(s) of the points that it holds after the distribution.

