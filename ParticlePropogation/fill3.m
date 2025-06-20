########################################################################

##

## Copyright (C) 2021 The Octave Project Developers

##

## See the file COPYRIGHT.md in the top-level directory of this

## distribution or <https://octave.org/copyright/>.

##

## This file is part of Octave.

##

## Octave is free software: you can redistribute it and/or modify it

## under the terms of the GNU General Public License as published by

## the Free Software Foundation, either version 3 of the License, or

## (at your option) any later version.

##

## Octave is distributed in the hope that it will be useful, but

## WITHOUT ANY WARRANTY; without even the implied warranty of

## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

## GNU General Public License for more details.

##

## You should have received a copy of the GNU General Public License

## along with Octave; see the file COPYING.  If not, see

## <https://www.gnu.org/licenses/>.

##

########################################################################


## -*- texinfo -*-

## @deftypefn  {} {} fill3 (@var{x}, @var{y}, @var{z}, @var{c})

## @deftypefnx {} {} fill3 (@var{x1}, @var{y1}, @var{z1}, @var{c1}, @var{x2}, @var{y2}, @var{z2}, @var{c2})

## @deftypefnx {} {} fill3 (@dots{}, @var{prop}, @var{val})

## @deftypefnx {} {} fill3 (@var{hax}, @dots{})

## @deftypefnx {} {@var{h} =} fill3 (@dots{})

## Create one or more filled 3-D polygons.

##

## The inputs @var{x}, @var{y}, and @var{z} are the coordinates of the polygon

## vertices.  If the inputs are matrices then the rows represent different

## vertices and each column produces a different polygon.  @code{fill3} will

## close any open polygons before plotting.

##

## The input @var{c} determines the color of the polygon.  The simplest form

## is a single color specification such as a @code{plot} format or an

## RGB-triple.  In this case the polygon(s) will have one unique color.  If

## @var{c} is a vector or matrix then the color data is first scaled using

## @code{caxis} and then indexed into the current colormap.  A row vector will

## color each polygon (a column from matrices @var{x}, @var{y}, and @var{z})

## with a single computed color.  A matrix @var{c} of the same size as @var{x},

## @var{y}, and @var{z} will compute the color of each vertex and then

## interpolate the face color between the vertices.

##

## Multiple property/value pairs for the underlying patch object may be

## specified, but they must appear in pairs.  The full list of properties is

## documented at @ref{Patch Properties}.

##

## If the first argument @var{hax} is an axes handle, then plot into this axes,

## rather than the current axes returned by @code{gca}.

##

## The optional return value @var{h} is a vector of graphics handles to

## the created patch objects.

##

## Example: oblique red rectangle

##

## @example

## @group

## vertices = [0 0 0

##             1 1 0

##             1 1 1

##             0 0 1];

## fill3 (vertices(:,1), vertices(:,2), vertices(:,3), "r");

## axis ([-0.5 1.5, -0.5 1.5, -0.5 1.5]);

## axis ("equal");

## grid ("on");

## view (-80, 25);

## @end group

## @end example

##

## @seealso{patch, fill, caxis, colormap}

## @end deftypefn


function h = fill3 (varargin)


  [hax, varargin] = __plt_get_axis_arg__ ("fill3", varargin{:});


  hlist = [];

  iargs = __find_patches__ (varargin{:});


  opts = {};

  if (numel (varargin) > iargs(end) + 3)

    opts = varargin(iargs(end)+4 : end);

  endif


  if (! all (cellfun (@(x) iscolorspec (x), varargin(iargs + 3))))

    print_usage ();

  endif


  oldfig = [];

  if (! isempty (hax))

    oldfig = get (0, "currentfigure");

  endif

  unwind_protect

    hax = newplot (hax);

    old_nxtplt = get (hax, "nextplot");

    if (! ishold ())

      set (hax, "box", "on");

    endif

    unwind_protect

      set (hax, "nextplot", "add");


      for i = 1 : length (iargs)

        x = varargin{iargs(i)};

        y = varargin{iargs(i) + 1};

        z = varargin{iargs(i) + 2};

        cdata = varargin{iargs(i) + 3};


        if (! size_equal (x, y, z))

          if (rows (x) != rows (y) || rows (x) != rows (z))

            error ("fill3: X, Y, and Z must have same number of rows");

          endif


          num_cols = max ([columns(x), columns(y), columns(z)]);

          if (iscolumn (x))

            x = repmat (x, [1, num_cols]);

          end

          if (iscolumn (y))

            y = repmat (y, [1, num_cols]);

          end

          if (iscolumn (z))

            z = repmat (z, [1, num_cols]);

          end

        endif


        if (isrow (x))

          x = x(:);

        endif

        if (isrow (y))

          y = y(:);

        endif

        if (isrow (z))

          z = z(:);

        endif


        if (ischar (cdata) || isequal (size (cdata), [1, 3]))

          one_color = true;

        else

          one_color = false;

        endif


        ## For Matlab compatibility, replicate cdata to match size of data

        if (! one_color && iscolumn (cdata))

          sz = size (x);

          if (all (sz > 1))

            cdata = repmat (cdata, [1, sz(2)]);

          endif

        endif


        ## For Matlab compatibility, return 1 patch object for each column

        for j = 1 : columns (x)

          if (one_color)

            [htmp, err] = __patch__ (hax, x(:,j), y(:,j), z(:,j), ...

                                     cdata, opts{:});

          else

            [htmp, err] = __patch__ (hax, x(:,j), y(:,j), z(:,j), ...

                                     cdata(:,j), opts{:});

          endif

          if (err)

            print_usage ();

          endif

          hlist(end+1, 1) = htmp;

        endfor

      endfor


      view (hax, 3);


    unwind_protect_cleanup

      if (strcmp (old_nxtplt, "replace"))

        set (hax, "nextplot", old_nxtplt);

      endif

    end_unwind_protect


  unwind_protect_cleanup

    if (! isempty (oldfig))

      set (0, "currentfigure", oldfig);

    endif

  end_unwind_protect


  if (nargout > 0)

    h = hlist;

  endif


endfunction


function iargs = __find_patches__ (varargin)

  iargs = 1:4:nargin;

  optidx = find (! cellfun (@isnumeric, varargin(iargs)), 1);

  iargs(optidx:end) = [];

endfunction


function retval = iscolorspec (arg)


  retval = false;

  if (ischar (arg))

    persistent colors = {"y", "yellow", "r", "red", "m", "magenta", ...

                         "c", "cyan", "g", "green", "b", "blue", ...

                         "w", "white", "k", "black"};

    if (any (strcmpi (arg, colors)))

      retval = true;

    endif

  elseif (isnumeric (arg))

    ## Assume any numeric argument is correctly formatted cdata.

    ## Let patch worry about the multiple different input formats.

    retval = true;

  endif


endfunction



%!demo

%! clf;

%! t1 = (1/16:1/8:1) * 2*pi;

%! t2 = ((1/16:1/8:1) + 1/32) * 2*pi;

%! x1 = sin (t1) - 0.8;

%! y1 = cos (t1);

%! z1 = sin (t1);

%! x2 = sin (t2) + 0.8;

%! y2 = cos (t2);

%! z2 = sin (t2);

%! h = fill3 (x1,y1,z1,"r", x2,y2,z2,"g");

%! title ({"fill3() function"; "cdata specified with string"});

%! grid ("on");


%!demo

%! clf;

%! t1 = (1/16:1/8:1) * 2*pi;

%! t2 = ((1/16:1/8:1) + 1/32) * 2*pi;

%! x1 = sin (t1) - 0.8;

%! y1 = cos (t1);

%! z1 = sin (t1);

%! x2 = sin (t2) + 0.8;

%! y2 = cos (t2);

%! z2 = sin (t2);

%! h = fill3 (x1,y1,z1,1, x2,y2,z2,2);

%! title ({"fill3() function"; 'cdata = row vector produces FaceColor = "flat"'});

%! grid ("on");


%!demo

%! clf;

%! x = [0 0

%!      1 0.5

%!      1 0.5

%!      0 0];

%! y = [0 0

%!      0 0

%!      1 0.5

%!      1 0.5];

%! z = y;

%! z(:,2) += 1e-4;

%! c = [1 2 3 4]';

%! fill3 (x, y, z, [c c]);

%! title ({"fill3() function"; 'cdata = column vector produces FaceColor = "interp"'});

%! grid ("on");

