<TeXmacs|WinTeXmacs-1.0.5>

<style|book>

<\body>
  <doc-data|<doc-title|BurrTools>|<doc-title|>|<doc-author-data|<author-name|Andreas
  R�ver>|<author-email|andreas_roever@web.de>>|<doc-note|<with-TeXmacs-text>>>

  \;

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Table of
    contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-1><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|1<space|2spc>The
    Program> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-2><vspace|0.5fn>

    1.1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-3>

    1.2<space|2spc><with|font-shape|small-caps|PuzzleSolver3D> Users
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-4>

    1.3<space|2spc>New Users <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-5>

    <with|par-left|1.5fn|1.3.1<space|2spc>Main Window
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-6>>

    <with|par-left|1.5fn|1.3.2<space|2spc>Piece Tab
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-8>>

    <with|par-left|1.5fn|1.3.3<space|2spc>Result Tab
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-10>>

    <with|par-left|1.5fn|1.3.4<space|2spc>Solve Tab
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-11>>

    1.4<space|2spc>Future Plans <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-12>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|2<space|2spc>The
    Internals> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-13><vspace|0.5fn>

    2.1<space|2spc>The puzzle file format
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-14>

    2.2<space|2spc>The Library <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-15>

    <with|par-left|1.5fn|2.2.1<space|2spc>Class voxel
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-16>>

    <with|par-left|1.5fn|2.2.2<space|2spc>Class puzzle
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-17>>

    <with|par-left|1.5fn|2.2.3<space|2spc>Class assembler
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-18>>

    <with|par-left|1.5fn|2.2.4<space|2spc>Class disassember
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-19>>

    <with|par-left|1.5fn|2.2.5<space|2spc>Example
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-20>>

    2.3<space|2spc>The Algorithms <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-21>

    <with|par-left|1.5fn|2.3.1<space|2spc>Assembly
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-22>>

    <with|par-left|1.5fn|2.3.2<space|2spc>Disassembly
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-23>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Epilogue>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <pageref|auto-24><vspace|0.5fn>
  </table-of-contents>

  <chapter|The Program>

  <section|Introduction>

  What are <name|BurrTools>? On the first hand they are two programs that
  assemble and disassemble burr-type puzzles. On is a contains a graphical
  interface while the other one is just as text program. On the other hand it
  is a C++ library that may help with the search and design of new puzzles of
  this type.

  The first section of the document is meant for the people already confident
  with <name|PuzzleSolver3D>. The second section is for the people completely
  new to this kind of software. And finally some internals of the program are
  explained.

  But first a little bit to the history of this program. Why another program
  you might ask. There are 2 main and one not so important reasons:

  <\enumerate-numeric>
    <item>The available programs (<name|PuzzleSolver3D> and
    <name|BCPBox/Genda>) are not for <name|Linux>, which is my operating
    system of choice.

    <item>The available programs are programs and wont allow me to do more
    interesting things like burr growing in an automated way.

    <item>The programs are quite expensive.
  </enumerate-numeric>

  Anyway, I was not satisfied with the available software. Then the C't, a
  German computer magazine, started a competition to write a program that
  counts the number of solutions to a merchandising puzzle of them as fast as
  possible. My program wasn't the fastest but it was the starting point for
  the <name|BurrTools>.

  As there are many people out there that are a lot more creative than I am
  and that could use a program like this to design nice puzzles I decided to
  make it public and free<\footnote>
    Free as in free speech and as in free beer (see
    <verbatim|http://www.gnu.org>)
  </footnote>.

  So I added a GUI that can work on many operating systems, including Linux
  and <name|Windows>. This has the disadvantage that the GUI looks a bit
  different from what the normal <name|Windows>-user is used to, so stay calm
  if things look a bit unusual, they behave in fact quite similar to what a
  normal <name|Windows>-program behaves.

  All this work has taken over two years to reach the current state, I hope
  it was worth it and you have a lot of fun with the program.

  <section|<name|PuzzleSolver3D> Users>

  People that know <name|PuzzleSover3D> from Andr� van Kamen should be able
  to quickly get used to the graphical user interface. For these people is
  this chapter. It just describes the important differences between
  <name|PuzzleSolver3D> and <name|BurrTools>. There are a few differences:

  <\enumerate-numeric>
    <item><name|BurrTools> doesn't handle holes automatically as
    <name|PuzzleSolver3D> does. This may at first sound like a disadvantage
    but in fact it is an advantage. <name|PuzzleSolver3D> normally teats all
    cubes of the target shape as cubes that might be filled but don't need to
    be. Except if you select ``Fill outer limits`` in the solve tab. Cubes
    that must be filled speed up the search process. The more there are, the
    faster the assembly will be as there are fewer possibilities.
    <name|BurrTools> requires for you to exactly specify which cubes in the
    result shape must be filled and which ones may be empty. This is done by
    either clicking with the left mouse key (must be filled cubes) or the
    right mouse button (may be filled cubes)

    <item><name|BurrTools> doesn't automatically find multiple identical
    pieces. You need to specify if there is one piece there more than once.
    If you just copy it the way you do in <name|PuzzleSolver3D> the program
    will find way too many solutions. For example for <name|Loveley18> it
    will find nearly 40'000'000 times as many solutions as there are in
    reality, so be careful.

    <item><name|BurrTools> allows you to define multiple problems in one
    file. So you can save all the SomaCube problems within one file

    <item><name|BurrTools> allows to put constrains on the placement of
    pieces. This is done by having an additional information attatched to
    each unit cube. I call this information color. You can then exactly
    define which colors in the pieces may go into which colored cubes inside
    the result. If just one of the conditions is not fullfilled the placement
    is not possible. To make thinks a bit easier there is a ``neutal'' called
    color that may go everywhere or accept every color.
  </enumerate-numeric>

  Then there are some other not so mayor things that differ:

  <\enumerate-numeric>
    <item>There are no limits to the number and size of pieces, like in
    <name|PuzzleSolver3D>. You can make the pieces as big as you want and are
    not limited to a size of 24

    <item>There is no limit to the number of placements to the pieces. It
    wont happen that the program comlains about too many placements. As long
    as you have enough memory the program will merrily continue working, even
    if it would take longer than the universe exists to complete the search.

    <item>You can not play with the puzzle. That is a useless feature (in my
    eyes)

    <item>You can not check the placements of the pieces

    <item>You don't need to clear solutions when you want to edit the puzzle
    again. The solutions get automatically removed as soon as you start
    another search or edit the puzzle.

    <item>There is no automatic zoom or animated rotation of the 3D view of
    the pieces or solution. Also the animation of the assembly process can
    not be animated (yet)

    <item>The solutions are not sorted by the number of moves required to
    disassemble the puzzle, so the hard solutions may be in the middle and
    hard to find.
  </enumerate-numeric>

  <section|New Users>

  This chapter is for users that are new to a program like this one. It will
  describe all aspects of the program in great detail.

  <subsection|Main Window>

  <\float|float|tbh>
    <big-figure||The Main Window<label|MainWindowImage>>
  </float>The main window (see Figure <reference|MainWindowImage>) contains a
  tool bar at the left side and the rest of the window is used for a 3
  dimensional view. By dragging with the mouse you can rotate the object
  visible inside the 3D view and by moving the slider at the very right edge
  you can zoom the object. The edge between the toolbar and the 3D view can
  be moved by dragging it with the mouse. Move the mouse to this edge. Once
  you see a left right array you can start to drag.

  The menu is relatively simple, just the usual ''New Puzzle'', \R''Load
  Puzzle'', \R''Save'' and \R''Save as'' entries.

  The toolbar on the left contain 3 main tabs. Each one of these tabs is
  explained in the following sections.

  <subsection|Piece Tab>

  <\float|float|tbh>
    <big-figure||The tools<label|ToolsImage>>
  </float>This tab (See Figure <reference|ToolsImage> left column) defines
  the colors and the shapes. <name|BurrTools> differentiates between shapes
  and pieces. There can be several pieces with the same shape. You as the
  user are responsible for not defining the same shape multiple times. If you
  do you will get multiple solutions.

  So let's start with the tab. At the top the color definition is placed
  ``Add'' and ``Remove'' add a new color or remove the selected color. Below
  is the list with the currently defined colors. There is always the
  ``Neutral'' color. This color can not be removed.

  Below the color definitions is the shape list followd by the shape editior.
  New adds another shape. \RDelete removes the current shape. \RCopy adds a
  new shape that is identical to the currently selected one.

  Below these 3 buttons is the list with the shapes. You can activate and
  edit a shape by clicking on it.

  Below the piece selector is another tab with some tools. These tools allow
  you to change the size of the space that is available for the piece
  definition. <name|BurrTools> uses the colours red green and blue to show
  the 3 axes. All tools that act on one of the axes are coloured accordingly.
  The 3D view contains three coloured lines that also show the axes. The next
  tab on the tool tab (see figure <reference|ToolsImage> right column bottom
  2 parts) contains buttons that allow you to rotate and shift and flip the
  shape definition inside its space. Finally the button minimise makes the
  definition space as small as possible. This is useful to make the files on
  disk a little bit smaller or to have a more centred view of the shape
  inside the 3D view. The Button ``Make inside variable'' sets all the unit
  cubes that are completely surrounded by cubes to variable cubes. If you
  apply this to your result shape the behaviour will be as in
  <name|PuzzleSolver3D>.

  The last item at the very bottom of this tab is the editor to change the
  shape. It contains a slider at the left that selects the z-plane. You edit
  by clicking into the squares shown. The left mouse button adds normal
  cubes. The right mouse button adds variable cubes. These variable cubes are
  only useful for result shapes. Clicking on an already filled cube either
  removes it or replaces it by the other cube type. Each cube changed will
  get the currently selected color (above in the color section) attatched to
  it. If this color it not the neutral color it will be visible as a small
  square in the upper left corner of the squares in the editor. The neutral
  color will not be visible anywhere.

  <twosuperior><subsection|Problem Tab>

  In this tab you assemble the problem(s) that you want to create with the
  shapes that you have defined. One problem consists of a result shape, some
  pieces that are supposed to assemble this shape and a list of color
  assignments.

  At the top is the list with your currenlty defined problems and some
  buttons to chreate new problems, delete them or copy an existing one.

  Below is the list with the problems. You can select a problem by clicking
  at it in your list. This problem is then edited.

  Below the problem list is the color assignment editor. You can see a list
  with the colors you have defined in the first tab, followed by some buttons
  folloed by the assignment list. Now this allows for every complex
  definitions, even thou most of the time very simple 1 to 1 assignments will
  be used.

  In the color assignment list you see the piece colors on the left and the
  result colors on the right connected by some arrows. The list is either
  sorted by the piece colors showing you then the piece color x can go into
  the following result colors, or the list is sorted by the result color
  showing you that piece colors x, y, ..., z can go into result color a. The
  2 Srt buttons switch between these views.

  The 2 buttons in the middle with the arrows add or remove a color from the
  currently selected entry in the color assignment list. The color ir always
  added on the side where the multiple colors are.

  Finally at the bottom is the shape assignment part. Here you select which
  shape is supposed to be assembled out of what other shapes. 2 List contain
  the shapes defined from the shapes tab and a list of the pieces involved in
  this problem. Above the shapelist is the result.

  Between the shape and the piece list are 2 buttons that either add another
  one of the selected pieces or remove one of it. In the lower list you can
  see how many of each piece are used for the problem.

  While editing problems you can see all the involved pieces in the 3D view.
  In the upper left coener the result is visible. It is drawn in double size.
  The upper right corner shows the currently selected piece and below these 2
  pieces are all the pieces visible that are included in the problem.

  <subsection|Solve Tab>

  Inside this tab (figure <reference|ToolsImage> 3rd column) you can start
  and stop the solving of a puzzle. But before you start there are a few
  options that may be useful. First the \Rdisassemble`` switch. If it is set
  the program tries to disassemble the assemblies found. Only the assemblies
  that disassemble are added to the list of solutions. Puzzles like the
  pentominos don't need this option set. It would only slow down the
  computations.

  Then there is the ''just count'' switch. Press this if you are not
  interested in the solutions themselves but only in the number of them.

  Finally the \Rreduce`` switch. This does some extra work before starting
  the real assembler process to find placements of pieces that are apparently
  not possible. I have not found a puzzle where this additional work at the
  start speeds up the search enough to regain the lost time. But sometimes
  this takes a lot of time and if you only want to see if there are solutions
  at all you can skip this part.

  The 3 buttons that start the search, continue a search that was stopped and
  stop a search. The program tries to solv the problem that is selected in
  the top list.

  Once a solution has been found it will be displayed in the 3D view. You can
  rotate it as usual.

  On the very bottom there is also a list of all the pieces visible. On this
  list the pieces that occur multiple times have an additional number
  appended to the number of the shape they belong to. Also the colour is
  slightly modified for each piece. If you click on the field for one piece
  you can toggle between normal view, wire-frame view and invisible mode for
  each piece. This is helpful if you want to shift pieces into a box. The box
  would hide most of the action that is going on inside it.

  <section|Future Plans>

  So, what are my future plans? There are a lot of things that are still
  missing from the current program. A list of things that might be
  interesting to implement are the following things:

  <\itemize-dot>
    <item>Include the burr growing into the GUI. This might not happen any
    time soon but it might be useful.

    <item>Include a possibility to generate images with solutions. Currently
    you need the program to animate for you.

    <item>Add some special algorithms that are faster for certain kind of
    puzzles. The current algorithms is quite good for nearly all puzzles, but
    it's not the fastest.

    <item>The disassembly could use some speedups. Here <name|PuzzleSolver3D>
    is sometimes much faster (but also sometimes much slower)
  </itemize-dot>

  I would be very happy to get contributions from other people. After all
  there are quite a few people out there that have their own puzzle solving
  programs, maybe they have some nice additions. There is one important thing
  to keep in mind: the additions have to run on <name|Linux>. So you can not
  use any proprietary library that is not available for <name|Linux>.

  <chapter|The Internals>

  <section|The puzzle file format>

  For those people that want to do things that the GUI is not supporting the
  exact file format of \ the files used by the GUI and the library may be of
  interest.

  The foamrt is actually a gzip compressed XML-File. The program can read
  both, compressed and uncompressed files transparently so you don't need to
  zip them bevore loading into the program.

  I wont describe all the elements of the XML-File, it's easier if you edit
  something similar inthe GUI and lokk in which way the program saves these
  information.

  <section|The Library>

  The library is available for all people who want to do an analysis that
  would be too much work to do by hand with the GUI. A bit of C++ programming
  experience is necessary to handle the task.

  There are 4 important classes in the library. The class <verbatim|voxel_c>
  handles a 3 dimensional array. Each position inside the array corresponds
  with one cube inside the piece. The class <verbatim|puzzle_c> is
  responsible for the whole puzzle containing a set of pieces and a solution.
  The classes <verbatim|assembler_x_c> and <verbatim|disassembler_x_c> (where
  <verbatim|x> is a number which may be available to select different
  algorithms that do the same task) are responsible to find assemblies and to
  disassemble the found assemblies. The important aspects of these classes
  will be explained in the next sections.

  <subsection|Class voxel>

  This class contains function so organize, modify, transform 3-dimensional
  arrays of cubes. There are 2 types of voxel arrays. One contains the
  definition for a single piece or a solution. The other contains an
  assembly. In the first type the entries inside the array can contain 3
  possible values:

  <\enumerate-numeric>
    <item><verbatim|VX_EMPTY>: the corresponding cube is empty

    <item><verbatim|VX_FILLED>: the corresponding cube is filled

    <item><verbatim|VX_VARIABLE>: the corresponding cube is unknown. These
    are the <em|may be filled> cubes inside solution shapes and growth area
    for the burrgrower.
  </enumerate-numeric>

  The class provides a set of functions to rotat, translate, mirror, resize
  and minimize the shape. The <verbatim|transform> function allows to
  generate all possible rotations <emdash> also including mirroring, if
  whished. The function <verbatim|selfSymmetries> calculates which of these
  transformations result in the same shape. <verbatim|Connected> finds out
  the all the cubes in the shape are connected in one big piece (neither the
  assembler nor the disassembler requests that this is the case).

  If all this is not enough then there are functions that return the value of
  the different cubes inside the shape and also to set the value of the
  cubes. These functions exists in different versions. One requires the x, y
  and z coordinate of the cube requestet. The other just takes one number.
  For this function all the cubes are in one long row. This function is
  efficient to use if all cubes are traversed and an action is done that is
  independent of the exact position of this cube inside the shape. Finally
  there is a set of get functions that also work with coordinated outside the
  box of the shape. These function always return <verbatim|VX_EMPTY> for
  cubes outside the bounding box.

  <subsection|Class puzzle>

  This class contains all the shapes that define a puzzle. It contains a list
  of shapes for the puzzle pieces. Assigned to each piece shape is the number
  of pieces that are available for this shape. The class also holds the shape
  the pieces are supposed to be assembled in.

  <subsection|Class assembler>

  As already explained this class tries to find assemblies for a puzzle. It
  uses the dancing link algorithm explained later.

  The caller is informed about found solution via a class that the caller has
  to provide. This class contains a function. This function is called for
  each found assembly with the found solution as parameter.

  The caller can then do whatever he pleases. He can just cound the number of
  solutions. By just increasing a counter. He can save the found solutions.
  He can analyse, if the found solution is disassembable. If the caller wants
  to keep the found solution he has to make of copy of the provided object,
  because it is reused by the assembler class for each found solution.

  <subsection|Class disassember>

  The disassembler is currently available in various versions. All these
  version use a similar algorithm. They are just various states of
  optimisation for more and more complex puzzles.\ 

  <subsection|Example>

  A very simple example can be found within the source code of the project.
  Check the burrTxt sources. They just check a few command line options, load
  the puzzle and then solve it, no fuzz with user interface, multi threaded
  application, ...

  <section|The Algorithms>

  There are only two algorithms of interest inside this program. One is the
  assembly algorithm. This one is based on the ``Dancing Link'' algorithm
  from D.E.Knuth. I needed to update thi algorithm in 2 ways:

  <\enumerate-numeric>
    <item>We require cubes that <em|may be filled> as well as cubes that
    <em|must be filled>. The original algorithm only provides the 2nd type of
    cubes.

    <item>We need to do something about multiple identical pieces. The
    original algorithm will find <with|mode|math|<big|prod><rsub|s\<in\>shapes>num(s)!>
    as many solutions as there really are.
  </enumerate-numeric>

  The 2nd interesting algorithm is the disassembler. This is mainly a bredth
  first tree search over all possible placements of the pieces.

  <subsection|Assembly>

  I will describe the only the basics for the origianl dancing link
  algorithm. For further information read the document available on Mr.
  Knuths webpage (<verbatim|http://www-cs-faculty.stanford.edu/~knuth/musings.html>).

  The algorithm represents the puzzle as a matrix. In this matrix the first
  columns represents the pieces and the last columns represent one cube of
  the result each.

  Each line of the matrix corresponds to one possible placement of one piece
  inside the result. The column of the piece and the columns the represent
  the places inside the colution that the piece occupies with the placement
  are 1 inside the matrix. All the other cells are 0.

  The search itself runs on this matrix. We search for a set of lines that
  all the lines contain exactly one 1 in each column. This means that each
  piece must be used and each cube in the result must be filled.

  The algorithm does 2 operations on the matrix:

  <\enumerate-numeric>
    <item>Cover column n and uncover column n. This means that the column is
    removed from the matrix and no longer taken into accound for the search.
    When a column is removed all the rows that contain a 1 in this column
    will also be removed.

    <item>Cover and uncover row n. This means that we select this row for the
    set of rows that we search. The row covering also removes and reincludes
    all the columns that contain a 1 in this row. On these columns operation
    1 is performed.
  </enumerate-numeric>

  The 2nd operation can be interpreted as. Taking one piece and putting it
  inside the result at one possible place. This results in the fact that a
  few cubes of the result don't need to be observed any longer and all
  placements of all other pieces that collide with this placement don't need
  to be checked further.

  The cover and uncover operations are the inversion of one another. If we
  first cover something and then uncover it again the matrix is in exactly
  the same state.

  The algorithm is now recoursively trying all possibilities. It selects one
  column and then tries covers all rows that contain a 1 in this columns and
  then calls itself.

  It finished when there are either no more columns left. Then we have found
  a solution or there is one column with no rows. Then we have found a dead
  end and backtrack.

  This algorithm is per se not dependent on square cubes it is not dependent
  on any shape. You only need to transfere your puzzle into the matrix. Even
  William Waites puzzles should be possible. But as the square and cubes are
  most common I have for now only implemented this transformation.

  Now to the changes that I have done to this basic algorithm. There is first
  the matter with the 2 types of cubes. This is easily solved by removing the
  columns of the cubes that <em|may be filled> from the list of columns that
  need to be covered. They are still in the matrix, they just don't <em|need>
  to be covered to find a solution.

  The 2nd problem was much harder. How handle multiple identical pieces? The
  solution that I finally implemented is to enforce an order. All pieces get
  a number and all the placements get a number. If we now have 2 identical
  pieces <with|mode|math|a> and <with|mode|math|b> with
  <with|mode|math|a\<less\>b> I force that the placement of
  <with|mode|math|a>, <with|mode|math|p(a)> is also smaller than the
  placement of <with|mode|math|b> so <with|mode|math|p(a)\<less\>p(b)>. This
  is done by always placing all identical pieces in one go. The moment the
  algorithm decides to place one of the pieces that occure multiple times it
  will also place all the others and always check that these have larger
  placement numbers.

  <subsection|Disassembly>

  The disassembly algorithm is a breadth first tree search. In this tree
  every node represents one possible relative position of the pieces. To find
  out what can be moved in this node the algorithm Bill Cuttler used for his
  6 Piece Burr analysis is used. Some optimisations were done that are
  (hopefully) speeding up the whole process.

  <section|Adding to the Library>

  There is currently one useful thing besides the normal improvements that
  might be added to the library: Other puzzle types. The assembly algorithms
  is so abstract that it can cope with many different types of assembly
  puzzles, as long as they have some kind of pattern. Currently the assembler
  only supports puzzles made out of cubes but there is nothing that prevents
  solving puzzle where the base unit is a hexagon. Of course the disassembler
  can not do work with this kind of puzzles.

  To add other geometries the assembler is split into 2 parts. The dancing
  link algorithm and the algorithm that prepares the matrix for the dancing
  link algorithm. This preparation part is called the frontend.
</body>

<\initial>
  <\collection>
    <associate|font|roman>
    <associate|language|british>
    <associate|page-breaking|optimal>
    <associate|page-medium|papyrus>
    <associate|page-orientation|portrait>
    <associate|par-columns|1>
    <associate|par-hyphen|professional>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate||<tuple|1.1|5>>
    <associate|MainWindowImage|<tuple|1.1|6>>
    <associate|ToolsImage|<tuple|1.2|6>>
    <associate|auto-1|<tuple|<uninit>|4>>
    <associate|auto-10|<tuple|1.3.3|7>>
    <associate|auto-11|<tuple|1.3.4|7>>
    <associate|auto-12|<tuple|1.4|8>>
    <associate|auto-13|<tuple|2|9>>
    <associate|auto-14|<tuple|2.1|9>>
    <associate|auto-15|<tuple|2.2|9>>
    <associate|auto-16|<tuple|2.2.1|9>>
    <associate|auto-17|<tuple|2.2.2|9>>
    <associate|auto-18|<tuple|2.2.3|10>>
    <associate|auto-19|<tuple|2.2.4|10>>
    <associate|auto-2|<tuple|1|5>>
    <associate|auto-20|<tuple|2.2.5|10>>
    <associate|auto-21|<tuple|2.3|10>>
    <associate|auto-22|<tuple|2.3.1|10>>
    <associate|auto-23|<tuple|2.3.2|11>>
    <associate|auto-24|<tuple|2.4|12>>
    <associate|auto-3|<tuple|1.1|5>>
    <associate|auto-4|<tuple|1.2|5>>
    <associate|auto-5|<tuple|1.3|6>>
    <associate|auto-6|<tuple|1.3.1|6>>
    <associate|auto-7|<tuple|1.1|6>>
    <associate|auto-8|<tuple|1.3.2|6>>
    <associate|auto-9|<tuple|1.2|6>>
    <associate|footnote-1|<tuple|1|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|The Main Window<label|MainWindowImage>|<pageref|auto-7>>

      <tuple|normal|The tools<label|ToolsImage>|<pageref|auto-9>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Table
      of contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>The
      Program> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      1.1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>

      1.2<space|2spc><with|font-shape|<quote|small-caps>|PuzzleSolver3D>
      Users <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>

      1.3<space|2spc>New Users <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>

      <with|par-left|<quote|1.5fn>|1.3.1<space|2spc>Main Window
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1.5fn>|1.3.2<space|2spc>Piece Tab
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1.5fn>|1.3.3<space|2spc>Problem Tab
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1.5fn>|1.3.4<space|2spc>Solve Tab
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      1.4<space|2spc>Future Plans <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>The
      Internals> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      2.1<space|2spc>The puzzle file format
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>

      2.2<space|2spc>The Library <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>

      <with|par-left|<quote|1.5fn>|2.2.1<space|2spc>Class voxel
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1.5fn>|2.2.2<space|2spc>Class puzzle
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|1.5fn>|2.2.3<space|2spc>Class assembler
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|1.5fn>|2.2.4<space|2spc>Class disassember
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|1.5fn>|2.2.5<space|2spc>Example
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      2.3<space|2spc>The Algorithms <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>

      <with|par-left|<quote|1.5fn>|2.3.1<space|2spc>Assembly
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|1.5fn>|2.3.2<space|2spc>Disassembly
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      2.4<space|2spc>Adding to the Library
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>
    </associate>
  </collection>
</auxiliary>