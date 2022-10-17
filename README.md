# Image-Processing-In-Microscopy
Course material for "Image Processing in Microscopy" at the Friedrich Schiller University in Jena in the winter term 2022/23

## Organisation
See [Moodle](https://moodle.uni-jena.de/course/view.php?id=32267) for more details about the course itself.

### Course Registration
Officially, over Friedolin. Remember to register for the lecture, the exercises and then at some point for the exam.

### Homework
<details>
    <summary>Details</summary>
    * Homework 01
        * Empty (HW01.jl)[HW01.jl]
</details>

### Q&A Seminars
<details>
    <summary>Details</summary>

    
</details>


## Code
To download the files, we recommend `git`:
```
git clone git@github.com:bionanoimaging/Image-Processing-In-Microscopy-2022-2023.git
```
Usually via a _git pull_ you can update the code. If anything goes wrong which you can't fix, clone it again to a new folder.


### Julia Installation
Download the recent version 1.8.1 on the [Julia Website](https://julialang.org/downloads/).

#### Editor
We recommend using [Visual Studio Code](https://www.julia-vscode.org/), especially install the Julia and git plugin for VSCode.

#### Documentation 
Also check out the [documentation](https://docs.julialang.org/en/v1/manual/performance-tips/). It is the best resource for julia because many other pages are outdated.

##### Cheatsheet
There is a [Cheatsheet](https://juliadocs.github.io/Julia-Cheat-Sheet/) available.

### Activate Environment
Open the downloaded source folder with VSCode.
We are going to use [Pluto.jl](https://github.com/fonsp/Pluto.jl) for the notebooks.
Open the Julia REPL:
Try to type:
```julia
julia> using Pluto # type y to install Pluto
 │ Package Pluto not found, but a package named Pluto is available from a registry. 
 │ Install package?
 │   (@v1.8) pkg> add Pluto 
 └ (y/n/o) [y]: y


julia> Pluto.run()

Opening http://localhost:1235/?secret=sdCsckRR in your default browser... ~ have fun!

Press Ctrl+C in this terminal to stop Pluto
```

A browser should open from where you can try to open a notebook.
