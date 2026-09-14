# Decision Making Models

## The Key Concept

Most “decision-making” models are **evidence accumulation models**: they
assume that evidence for a decision or choice is accumulated over time
until it reaches a threshold. The time it takes to reach the threshold
is the observed response time (RT). Many “variants” of such models
exist, differing in the assumptions they make about the accumulation
process.

Wald Model

`drift` drift rate

``` js
// Both ends of the range are set by the frame, and they move with the boundary
// count. With one boundary every accumulation has to reach it inside the
// figure, which puts the floor at 1.5: the mean path then lands at 0.93 s in
// the slowest corner (boundary 0.8, ndt 0.4). With two, the far boundary takes
// a share 1 / (1 + exp(drift * boundary)) of the responses, which at a
// drift of 3 is one response in twenty - too few to be worth a curve. So the
// range moves down to where both distributions can be seen (0.18 of them at
// the far boundary at the default) and stops at 1.2, which is where the mean
// path still lands inside the frame in that same slowest corner.
viewof drift = Inputs.range(nbounds === 1 ? [1.5, 6] : [1.2, 4],
                            {value: nbounds === 1 ? 3 : 1.5, step: 0.1,
                             width: 190})
```

`boundary` decision threshold

``` js
// The range doubles with the boundary count, because what `boundary` measures
// does (see `edge`), while the figure stays exactly as it was: a boundary of 1
// with two of them draws the same two lines as 0.5 did as a threshold.
viewof boundary = Inputs.range(nbounds === 1 ? [0.3, 0.8] : [0.6, 1.6],
                               {value: nbounds === 1 ? 0.5 : 1,
                                step: nbounds === 1 ? 0.05 : 0.1, width: 190})
```

`bias` start point

``` js
// Only the two-boundary model has a start point to place, so `model` below
// shows and hides this one with the second boundary - it starts out hidden in
// the markup, so that it cannot be seen before the first cell runs, and it
// keeps its value while it is away. The range stops well short of either end:
// at `bias` 0.2 the start point is already only a fifth of the way up from the
// lower boundary, which is as lopsided as the figure can draw and still
// measure it.
viewof bias = Inputs.range([0.2, 0.8], {value: 0.5, step: 0.05, width: 190})
```

`ndt` non-decision time

``` js
viewof ndt = Inputs.range([0, 0.4], {value: 0.2, step: 0.01, width: 190})
```

``` js
// Both wrappers only decorate the figure and hand it back, so what follows
// reads as a plain `Plot.plot` call - dragging and animation are bolted on
// from the outside and nothing in the mark list knows about either.
dragHandles(fireTrials(Plot.plot({
  width: plotWidth,
  height: frame.height,
  marginLeft: cfg.mleft,
  marginRight: cfg.mright,
  marginTop: cfg.mtop,
  marginBottom: cfg.mbottom,
  // No ticks anywhere on time - the axis is the arrow drawn at the start
  // point below, and the one number worth reading is the mean RT.
  x: {domain: [0, cfg.tmax], axis: null},
  // Two numbers on the evidence axis, the start point and the boundary - three
  // once the boundary has a mirror image - with no tick marks: nothing else
  // here has a protruding tick, and the labels sit close enough to the frame
  // to read without one.
  // The label is rotated a quarter turn anticlockwise, which turns the arrow
  // it is written with at each end into one pointing down and one pointing up.
  // Plot's own `labelArrow` would add a third.
  y: {domain: [frame.ybot, frame.ytop], label: "← Evidence →",
      labelArrow: "none", labelAnchor: "center", ticks: frame.ticks,
      tickFormat: d3.format(".2f"), tickSize: 0, tickPadding: 7},
  marks: [
    // The RT distribution, sitting on the boundary it is the crossing time of.
    Plot.areaY(density.up, {x: "t", y1: "base", y2: "y", fill: model.colour,
                            fillOpacity: 0.15}),
    Plot.line(density.up, {x: "t", y: "y", stroke: model.colour,
                           strokeWidth: 2}),
    // The other boundary's distribution, hanging under it - empty when there
    // is only one boundary. The two are drawn to a single scale, so the lower
    // one is shorter by exactly the share of responses it takes. That share is
    // what a second boundary adds to the model, and giving each curve a height
    // of its own would throw it away.
    Plot.areaY(density.low, {x: "t", y1: "base", y2: "y", fill: model.colour,
                             fillOpacity: 0.15}),
    Plot.line(density.low, {x: "t", y: "y", stroke: model.colour,
                            strokeWidth: 2}),
    // Thirty trials, redrawn from scratch whenever a slider moves.
    Plot.line(trials.rows, {x: "t", y: "x", z: "id", stroke: model.colour,
                            strokeOpacity: 0.42, strokeWidth: 0.9, clip: true}),
    // Time runs along the start point, and rides up and down with it. Black,
    // because it is the one line here that is an axis rather than an
    // annotation.
    Plot.arrow([{x1: 0, y1: start, x2: cfg.tmax, y2: start}],
               {x1: "x1", y1: "y1", x2: "x2", y2: "y2", stroke: cfg.black,
                strokeWidth: 1.2, headLength: 7}),
    Plot.text([{t: cfg.tmax, y: start}],
              {x: "t", y: "y", text: ["Time"], dx: -12, dy: 13,
               textAnchor: "end", fill: cfg.black, fontSize: 11}),
    // Nothing accumulates before ndt. The double arrow lies along the time
    // arrow, over the flat stretch every path starts with, and is the whole
    // annotation: a dropped line up to the boundary only repeated what the
    // foot of the drift arrow already marks.
    Plot.arrow(ndtSpan, {x1: "x1", y1: "y1", x2: "x2", y2: "y2",
                         stroke: cfg.purple, strokeWidth: 1.8, headLength: 8}),
    Plot.text(ndtSpan.slice(0, 1),
              {x: "x1", y: "y1", text: ["ndt"], dy: -8, fill: cfg.purple,
               fontSize: 11, stroke: "white", strokeWidth: 3,
               paintOrder: "stroke"}),
    // How far up from the lower boundary the evidence starts, which is what
    // `bias` is. At the left edge, where nothing has happened yet - and empty
    // with one boundary, which has no start point to place.
    Plot.arrow(biasSpan, {x1: "x1", y1: "y1", x2: "x2", y2: "y2",
                          stroke: cfg.teal, strokeWidth: 1.8, headLength: 8}),
    Plot.text(biasSpan.slice(0, 1),
              {x: "x1", y: "y1", text: ["bias"], dx: 5, textAnchor: "start",
               fill: cfg.teal, fontSize: 11, stroke: "white", strokeWidth: 3,
               paintOrder: "stroke"}),
    // The boundary is the second axis of the figure, so it is drawn at the
    // weight of one rather than laid over the frame as a heavier line. Its
    // mirror image carries no label of its own: the tick states its value, and
    // the two lines are one parameter.
    Plot.ruleY(frame.rules, {stroke: cfg.orange, strokeWidth: 1.2}),
    Plot.text([{t: cfg.tmax, y: edge}],
              {x: "t", y: "y", text: ["boundary"], dy: -7, textAnchor: "end",
               fill: cfg.orange, fontSize: 11}),
    // The boundary is a distance, not a place on the clock, and which distance
    // depends on how many boundaries there are. One dimension line at the right
    // edge, under the label, spans it: from the start point up to the threshold
    // with one boundary, and from one boundary to the other with two.
    Plot.arrow(boundarySpan, {x1: "x1", y1: "y1", x2: "x2", y2: "y2",
                              stroke: cfg.orange, strokeWidth: 1.8,
                              headLength: 8}),
    // The mean drift path, from the start point: its slope is the drift rate,
    // and where it lands is the time a noiseless accumulation would get there.
    // The arc at its foot marks that slope as an angle.
    Plot.arrow([driftPath],
               {x1: "x1", y1: "y1", x2: "x2", y2: "y2", stroke: cfg.green,
                strokeWidth: 3, headLength: 11}),
    Plot.arrow(driftArc.data, {x1: "x1", y1: "y1", x2: "x2", y2: "y2",
                               stroke: cfg.green, strokeWidth: 1.8,
                               headLength: 8, bend: driftArc.bend}),
    // The white stroke is a halo: these labels sit over the path tangle. The
    // label sits in the wedge the arc marks out rather than partway up the
    // arrow, which at a steep drift would put it on top of the boundary.
    Plot.text([driftArc.label],
              {x: "x", y: "y", text: ["drift"], dx: 2, textAnchor: "start",
               fill: cfg.green, fontSize: 11, stroke: "white", strokeWidth: 3,
               paintOrder: "stroke"}),
    Plot.dot(trials.hits, {x: "t", y: "x", fill: cfg.orange, r: 3}),
    Plot.ruleX([ndt + meanDT],
               {y1: edge, y2: frame.ytop, stroke: cfg.grey,
                strokeDasharray: "4,3"}),
    Plot.text([{t: ndt + meanDT}],
              {x: "t", text: ["mean RT " + d3.format(".2f")(ndt + meanDT) + " s"],
               frameAnchor: "top", textAnchor: "start", dx: 4, dy: 8,
               fill: cfg.grey, fontSize: 11})
  ]
})))
```

Assumptions

`boundaries`

``` js
viewof nbounds = Inputs.radio([1, 2], {value: 1})
```

``` js
// tmax and the slider ranges are set against each other. Both boundary ranges
// put the lines at the same place - an `edge` of 0.3 to 0.8 - so the figure is
// sized against that rather than against the slider: the slowest corner of the
// parameter space (drift at the floor of its range, edge 0.8, ndt 0.4) has a
// mean RT of 0.94 s, so 1.2 s of axis holds every setting without the frame
// ever having to rescale - which it must not do, or raising ndt would stop
// looking like a shift. Same for the evidence axis: ytop clears the tallest
// density (edge 0.8 + dheight) and the frame is otherwise as tight as it can
// be, because the accumulation band is only as tall as `edge` and the smallest
// one has to stay readable - which is why neither range goes below it. ybot
// leaves room for the paths, which wander below zero on the way up; with two
// boundaries they cannot go past the second one and the frame turns symmetric
// instead - see `frame`, which owns everything the boundary count moves,
// heights included. driftcap is how far short of the right edge the drift
// arrow stops when the boundary is further off than the figure is wide - a
// tenth of a second, which leaves the arrowhead clear of the boundary label in
// that corner.
//
// The colours are read back out of the `:root` block at the top of this file
// rather than written again here, so the palette is stated once and the
// figure cannot fall out of step with the controls. A name that is not
// declared up there throws rather than coming back as an empty string, which
// would otherwise paint a mark in nothing at all and say nothing about why.
cfg = {
  const css = getComputedStyle(document.documentElement);
  const hex = (name) => {
    const v = css.getPropertyValue("--cogmod-" + name).trim();
    if (!v) throw new Error("--cogmod-" + name + " is not in the palette");
    return v;
  };
  return {tmax: 1.2, dt: 0.005, npaths: 30, ybot: -0.5, ytop: 1.4,
          dheight: 0.55, driftcap: 0.1, mleft: 46, mright: 18, mtop: 18,
          mbottom: 16, blue: hex("blue"), azure: hex("azure"),
          orange: hex("orange"), green: hex("green"), teal: hex("teal"),
          purple: hex("purple"), grey: hex("grey"), black: hex("black")};
}
```

``` js
// Where the boundary lines are drawn, either side of the middle of the figure.
// With one boundary that is `boundary` itself, the distance the evidence has
// to cover; with two, `boundary` is the whole gap between them, so each line
// sits at half of it. Everything geometric reads this rather than the slider -
// the lines, the ticks, the paths' stopping points, the drift arrow's target -
// while the two densities, which take the separation as such, read `boundary`.
edge = nbounds === 1 ? boundary : boundary / 2
```

``` js
// The start point, in evidence units, and so the height the time axis is drawn
// at. With one boundary it is the zero of the axis and stays there. With two,
// `bias` slides it between them - the same proportion `dcogmod_ddm()` takes,
// measured from the lower boundary, so 0.5 starts midway and 0.75 three
// quarters of the way up. Everything that begins where the evidence begins
// reads it: the time arrow, the ndt dimension line that now sits on that
// arrow, the foot of the drift arrow, and the paths.
start = nbounds === 1 ? 0 : edge * (2 * bias - 1)
```

``` js
// The evidence axis, which is the one thing the boundary count moves. Two
// boundaries make it symmetric: either curve can be the taller of the two -
// a start point placed low makes the lower boundary the likelier one - so the
// frame has to hold a full-height density on both sides, which is `ytop` at
// each end. The figure is given the extra height rather than the extra domain
// alone, so a unit of evidence stays the same number of pixels tall as it was
// (140 against 140) and nothing looks squashed by the second boundary - which
// is the whole reason the two heights are stated here, next to the domains
// they have to be kept in proportion with, rather than among the constants in
// `cfg`. The ticks name the two boundaries and the start point between them;
// the start point is the only one of the three that is not also drawn as a
// line, the time axis being that line.
frame = nbounds === 1
  ? {ybot: cfg.ybot, ytop: cfg.ytop, height: 300,
     ticks: [0, edge], rules: [edge]}
  : {ybot: -cfg.ytop, ytop: cfg.ytop, height: 426,
     ticks: [-edge, start, edge], rules: [-edge, edge]}
```

``` js
// Two boundaries make this a different model, and everything that says which
// one is on show lives here: the name the title takes, and the colour the
// density and the paths are drawn in - which the chosen segment of the toggle
// wears too, so that the switch says what it is about to do. The heading is
// written in the markdown above rather than returned from here, so that the
// column it heads keeps the spacing it has; this cell only fills in the name.
model = {
  const it = nbounds === 1 ? {name: "Wald Model", colour: cfg.blue}
                           : {name: "Drift Diffusion Model", colour: cfg.azure};
  const el = document.querySelector(".cogmod-title");
  if (el) el.textContent = it.name;
  // `bias` has nothing to say about a single boundary, so its slider leaves
  // with the second one.
  const biased = document.querySelector(".cogmod-par-bias");
  if (biased) biased.hidden = nbounds === 1;
  return it;
}
```

``` js
// `.cogmod-controls` takes 200px plus the 1.75rem flex gap; Inputs.range does
// not carry its own label here, so the column is exactly as wide as it says.
plotWidth = Math.max(320, Math.min(width, 900) - 232)
```

``` js
// Three of the marks are handles. The boundary line drags up and down, an
// invisible vertical line at ndt drags left and right, and the drift arrow
// pivots about its foot - each writing to its own slider, so a value still
// lives in exactly one place and the figure is only a second way in.
//
// Every change rebuilds the plot, so the element the pointer went down on is
// gone by the second frame of a drag. The move and release listeners therefore
// go on the window, and the geometry they work from is measured once at
// pointerdown - which stays true, because each rebuild is the same size in the
// same place. The figure is drawn in its own svg units, so a screen pixel is
// worth `box.height / rect.height` of them: the two agree at full size and
// part company when the column is narrow enough for Plot's max-width to kick
// in.
dragHandles = function(plot) {
  const GRAB = 7;                    // how near the pointer has to come, in px
  const slider = {
    drift: (viewof drift).querySelector('input[type="range"]'),
    boundary: (viewof boundary).querySelector('input[type="range"]'),
    bias: (viewof bias).querySelector('input[type="range"]'),
    ndt: (viewof ndt).querySelector('input[type="range"]')
  };
  const at = (key) => Number(slider[key].value);
  // Where the lines are and where the evidence starts, read off the sliders
  // rather than taken from `edge` and `start`, for the same reason every other
  // handle reads its own value live.
  const edgeAt = () => nbounds === 1 ? at("boundary") : at("boundary") / 2;
  const origin = () => nbounds === 1 ? 0 : edgeAt() * (2 * at("bias") - 1);
  const box = plot.viewBox.baseVal;
  const perPixel = (rect) => box.height / rect.height;
  // Where a value sits on screen, and the value at a point on screen: the same
  // map read each way, once per axis.
  const mapper = (scale, anchor) => {
    const [d0, d1] = scale.domain, [p0, p1] = scale.range;
    const k = (p1 - p0) / (d1 - d0);
    return {
      to: (v, rect) => anchor(rect) + (p0 + (v - d0) * k) / perPixel(rect),
      from: (p, rect) => d0 + ((p - anchor(rect)) * perPixel(rect) - p0) / k
    };
  };
  const X = mapper(plot.scale("x"), (rect) => rect.left);
  const Y = mapper(plot.scale("y"), (rect) => rect.top);
  const toSegment = (p, a, b) => {
    const vx = b.x - a.x, vy = b.y - a.y, len = vx * vx + vy * vy;
    const t = len ? Math.max(0, Math.min(1, ((p.x - a.x) * vx +
                                             (p.y - a.y) * vy) / len)) : 0;
    return Math.hypot(p.x - (a.x + t * vx), p.y - (a.y + t * vy));
  };

  // `away` is the distance from the pointer to the handle; `read` is the value
  // the pointer is asking that parameter to take.
  const handle = {
    boundary: {
      cursor: "ns-resize",
      // Either line is the handle when there are two of them, and the lower
      // one asks for the same value the upper would - one parameter drawn
      // twice - which is what the distance from the start point reads off.
      away: (p, rect) => {
        const b = edgeAt();
        return Math.min(...(nbounds === 1 ? [b] : [b, -b])
          .map((y) => Math.abs(p.y - Y.to(y, rect))));
      },
      // Once there are two, the line dragged is one side of a gap that is
      // symmetric about the middle, so the value it asks for is twice its
      // distance from there - and a line dragged through the middle comes out
      // the other side rather than sticking at the bottom of the slider.
      read: (p, rect) => nbounds === 1 ? Y.from(p.y, rect)
                                       : 2 * Math.abs(Y.from(p.y, rect))
    },
    ndt: {
      cursor: "ew-resize",
      // Nothing is drawn at ndt any more, so the line to grab is the one the
      // double arrow ends on, and only over the stretch of the figure ndt has
      // anything to say about: the accumulation band, between the boundaries
      // the evidence is on its way to.
      away: (p, rect) =>
        p.y < Y.to(edgeAt(), rect) - GRAB ||
        p.y > Y.to(nbounds === 1 ? 0 : -edgeAt(), rect) + GRAB
          ? Infinity : Math.abs(p.x - X.to(at("ndt"), rect)),
      read: (p, rect) => X.from(p.x, rect)
    },
    drift: {
      cursor: "move",
      // The arrow turns about its foot, so the rate the pointer asks for is
      // the slope of the line from the foot out to it - evidence over time,
      // which is what a drift rate is.
      away: (p, rect) => {
        const o = origin();
        const tip = driftTip(at("ndt"), edgeAt(), o, at("drift"));
        return toSegment(p,
          {x: X.to(at("ndt"), rect), y: Y.to(o, rect)},
          {x: X.to(tip, rect),
           y: Y.to(o + at("drift") * (tip - at("ndt")), rect)});
      },
      read: (p, rect) => (Y.from(p.y, rect) - origin()) /
                         (X.from(p.x, rect) - at("ndt"))
    },
    bias: {
      cursor: "ns-resize",
      // The start line itself, which is the time axis: dragging it slides the
      // start point between the boundaries. Only to the right of ndt, where
      // the line is free - to the left of it the ndt dimension arrow lies on
      // top and its own handle stands there.
      away: (p, rect) => nbounds === 1 ||
                         p.x < X.to(at("ndt"), rect) + GRAB
        ? Infinity : Math.abs(p.y - Y.to(origin(), rect)),
      read: (p, rect) => (Y.from(p.y, rect) + edgeAt()) / (2 * edgeAt())
    }
  };

  // Lines this thin are too fine to aim at, so anything within GRAB counts,
  // and where two handles are both in reach the nearer one wins.
  const nearest = (event, rect) => {
    const p = {x: event.clientX, y: event.clientY};
    let pick = null, best = GRAB;
    for (const key in handle) {
      const away = handle[key].away(p, rect);
      if (away < best) { pick = key; best = away; }
    }
    return pick;
  };

  plot.addEventListener("pointermove", (event) => {
    const pick = nearest(event, plot.getBoundingClientRect());
    plot.style.cursor = pick ? handle[pick].cursor : "";
  });

  plot.addEventListener("pointerdown", (event) => {
    const rect = plot.getBoundingClientRect();
    const pick = nearest(event, rect);
    if (!pick) return;
    event.preventDefault();
    document.body.style.userSelect = "none";
    const input = slider[pick];
    const lo = Number(input.min), hi = Number(input.max);
    const step = Number(input.step);
    const move = (e) => {
      const asked = handle[pick].read({x: e.clientX, y: e.clientY}, rect);
      if (!Number.isFinite(asked)) return;      // straight above the pivot
      const v = Math.max(lo, Math.min(hi, Math.round(asked / step) * step));
      if (v === Number(input.value)) return;
      input.value = v;
      input.dispatchEvent(new Event("input", {bubbles: true}));
    };
    const stop = () => {
      window.removeEventListener("pointermove", move);
      window.removeEventListener("pointerup", stop);
      document.body.style.userSelect = "";
    };
    window.addEventListener("pointermove", move);
    window.addEventListener("pointerup", stop);
  });

  return plot;
}
```

``` js
// A fresh set of trials every seven seconds - long enough for the one before
// it to have finished arriving and stood complete for a beat. The figure is an
// explainer, so it runs on its own - but not for a reader who has asked their
// machine for less movement, who gets one static set and no timer at all, and
// not while the figure is off screen, because this sits near the top of a long
// article and a figure nobody is looking at has no business simulating
// anything. The ndt slider is watched rather than the plot: it is in the same
// row, and unlike the plot - and unlike the drift and boundary sliders, which
// the boundary count rebuilds with new ranges - it is never replaced out from
// under the observer.
volley = window.matchMedia("(prefers-reduced-motion: reduce)").matches ? null
  : Generators.observe((next) => {
      let n = 0, timer = null;
      next(n);
      const watch = new IntersectionObserver(([seen]) => {
        if (seen.isIntersecting && timer === null) {
          timer = setInterval(() => next(++n), 7000);
        } else if (!seen.isIntersecting && timer !== null) {
          clearInterval(timer);
          timer = null;
        }
      });
      watch.observe(viewof ndt);
      return () => {
        if (timer !== null) clearInterval(timer);
        watch.disconnect();
      };
    })
```

``` js
// True once per volley. The plot is rebuilt whenever a parameter changes too,
// and those rebuilds must not restart the sweep: the traces would be wiped off
// the figure for as long as a slider or a handle was moving. Cell state, not
// cell value - nothing depends on it, so it is created once and remembers.
newVolley = {
  let seen = null;
  return (n) => (n === seen ? false : (seen = n, true));
}
```

``` js
// The trials are fired one at a time rather than arriving as a block: each
// trace gets a clip of its own and they are let go a fifth of a second apart,
// so the reader watches evidence pile up thirty times over instead of once.
// Each trace draws at three times the model's own speed - at real time a
// single one takes a second, which is long enough that the set stops reading
// as one thing - and the dot a trace ends on rides that same trace's clip, so
// it lands as the trace reaches the boundary rather than when the set is done.
//
// A set takes 29 x 200 + 400 = 6.2 s to fill, and the timer below comes round
// at 7 s: the last trace holds for a moment, then the figure clears and starts
// over. Clearing is the cheap end of the two - keeping the traces and dropping
// the oldest would mean a rolling buffer, a launch time carried per trace, and
// a rebuild every time one was fired.
//
// Each rect is full width and *scaled* down to nothing, rather than drawn at
// zero width and grown: a browser that will not run the animation is then left
// with rects that cover everything and traces that are simply all there.
fireTrials = function(plot) {
  if (volley === null || !newVolley(volley)) return plot;
  const SPEED = 3;                       // times the model's own clock
  const STAGGER = 200;                   // ms between one trace and the next
  const svgns = "http://www.w3.org/2000/svg";
  const box = plot.viewBox.baseVal;
  const [x0, x1] = plot.scale("x").range;
  // The trials are the one line mark drawing more than a single path; the rest
  // of the figure is the frame, and stays.
  const held = [...plot.querySelectorAll('g[aria-label="line"]')]
    .find((g) => g.querySelectorAll("path").length > 1);
  if (!held) return plot;
  const clipOf = [...held.querySelectorAll("path")].map((path, i) => {
    const id = "cogmod-shot-" + volley + "-" + i;
    const clip = document.createElementNS(svgns, "clipPath");
    clip.setAttribute("id", id);
    const front = document.createElementNS(svgns, "rect");
    front.setAttribute("x", x0);
    front.setAttribute("y", 0);
    front.setAttribute("width", x1 - x0);
    front.setAttribute("height", box.height);
    front.style.transformBox = "view-box";
    front.style.transformOrigin = x0 + "px 0px";
    clip.appendChild(front);
    plot.insertBefore(clip, plot.firstChild);
    path.setAttribute("clip-path", "url(#" + id + ")");
    // `both` holds the first frame through the delay, so a trace waiting its
    // turn is clipped away rather than sitting there whole.
    front.animate([{transform: "scaleX(0)"}, {transform: "scaleX(1)"}],
                  {duration: cfg.tmax * 1000 / SPEED, delay: i * STAGGER,
                   easing: "linear", fill: "both"});
    return id;
  });
  // Plot draws the dots in the order the hits were recorded, and the paths one
  // per trial in trial order, which is what lets a dot find its own trace.
  const dots = plot.querySelector('g[aria-label="dot"]');
  if (dots) {
    [...dots.children].forEach((dot, i) => {
      const hit = trials.hits[i];
      if (hit && clipOf[hit.id]) {
        dot.setAttribute("clip-path", "url(#" + clipOf[hit.id] + ")");
      }
    });
  }
  return plot;
}
```

``` js
// The value bubble each slider shows while the pointer is on it. This runs
// once per input rather than on every move: it names the four `viewof`
// elements, not their values, so ojs does not re-run it when a slider moves -
// the listener takes over from there. Showing and hiding is left to CSS.
//
// The thumb is 13px wide and slides between the two ends of the track, so its
// centre travels the width of the input less one thumb - which is the offset
// the bubble has to follow to sit over it.
bubbles = {
  for (const form of [viewof drift, viewof boundary, viewof bias,
                      viewof ndt]) {
    // This runs again whenever the boundary count rebuilds one of these. The
    // ones it did not rebuild are the same elements as before, bubble and
    // listeners and all, so they are left alone: replacing the bubble there
    // would strand the old one with a live listener still writing to it.
    if (form.querySelector(".cogmod-bubble")) continue;
    const range = form.querySelector('input[type="range"]');
    const bubble = form.appendChild(document.createElement("span"));
    bubble.className = "cogmod-bubble";
    const place = () => {
      const lo = Number(range.min), hi = Number(range.max);
      const frac = (range.valueAsNumber - lo) / (hi - lo);
      bubble.textContent = range.value;
      bubble.style.left = (range.offsetLeft + 6.5 +
                           frac * (range.offsetWidth - 13)) + "px";
    };
    range.addEventListener("input", place);
    // A slider that was hidden when this ran has no width to measure, so the
    // bias one is placed again as the pointer arrives - which is the moment
    // before the bubble is shown.
    form.addEventListener("pointerenter", place);
    place();
  }
}
```

``` js
// Two arrows pointing outwards from the midpoint, which is where the label
// goes. Below about 40 ms of ndt the heads would be longer than the span they
// measure, so the whole annotation drops out rather than turn into a blob.
ndtSpan = ndt > 0.04
  ? [{x1: ndt / 2, y1: start, x2: 0, y2: start},
     {x1: ndt / 2, y1: start, x2: ndt, y2: start}]
  : []
```

``` js
// The same trick standing up at the left edge, measuring the start point from
// the lower boundary - which is `bias`, times the separation. Nothing else is
// drawn out here: the paths do not start until ndt, and the densities are flat
// against their boundaries this early. The dropout guard is the ndt one's,
// in evidence units: below about a sixth of one the two heads would be longer
// than the span between them.
biasSpan = {
  const bx = 0.025;
  const z = start + edge;              // up from the lower boundary
  return nbounds === 1 || z < 0.16 ? []
    : [{x1: bx, y1: start - z / 2, x2: bx, y2: -edge},
       {x1: bx, y1: start - z / 2, x2: bx, y2: start}];
}
```

``` js
// The same two-arrows-from-the-midpoint trick standing up, at the right edge
// of the frame where the paths have all finished. One arrow for one parameter:
// it spans whatever stretch the slider sets (see `edge`). Neither span moves
// with `bias` - the first is measured at a start point that cannot move, and
// the second is between the boundaries, which the start point slides between
// without changing. No dropout guard, because the shortest boundary the slider
// allows is still 40-odd pixels tall. bx is a hair in from tmax so the strokes
// clear the edge.
boundarySpan = {
  const bx = cfg.tmax - 0.02;
  const [lo, hi] = nbounds === 1 ? [0, edge] : [-edge, edge];
  const mid = (lo + hi) / 2;
  return [{x1: bx, y1: mid, x2: bx, y2: lo},
          {x1: bx, y1: mid, x2: bx, y2: hi}];
}
```

``` js
// Where the mean drift arrow ends: at the upper boundary, or `driftcap` short
// of the right edge when a slow drift starting low down would need more of the
// clock than the figure has to get there. What the arrow draws is the slope,
// and the slope is the same either way. Written here rather than at either of
// its two callers - `driftPath` draws the arrow and the drift handle measures
// the pointer against it, off live slider values - so that the cap cannot end
// up meaning one thing to the figure and another to the drag.
driftTip = function(ndt, edge, start, drift) {
  return Math.min(ndt + (edge - start) / drift, cfg.tmax - cfg.driftcap);
}
```

``` js
// The mean drift path: out of the start point at the drift rate, as far as
// `driftTip` allows.
driftPath = {
  const x2 = driftTip(ndt, edge, start, drift);
  return {x1: ndt, y1: start, x2: x2, y2: start + drift * (x2 - ndt)};
}
```

``` js
// The arc at the foot of the drift arrow, marking the angle it makes with
// time. Both ends have to be the same distance from the corner *on screen*,
// and the two axes carry different units, so this measures the frame in
// pixels and converts back. The radius shrinks with the arrow so that a steep,
// short arrow does not end up with an arc longer than itself.
//
// `bend` is the angle between the chord and the tangent at its ends, which for
// a circular arc is half the angle it subtends - so passing half the drift
// angle is what makes this an arc centred on the corner rather than a line
// cutting it off. It has to be computed: at 40 degrees, Plot's own default
// bend of 22.5 leaves a sagitta of two pixels and the curve reads as straight.
driftArc = {
  const pxPerSec = (plotWidth - cfg.mleft - cfg.mright) / cfg.tmax;
  const pxPerEv = (frame.height - cfg.mtop - cfg.mbottom) /
                  (frame.ytop - frame.ybot);
  const dx = (driftPath.x2 - driftPath.x1) * pxPerSec;
  const dy = (driftPath.y2 - driftPath.y1) * pxPerEv;
  const len = Math.hypot(dx, dy);
  const r = Math.min(34, 0.45 * len);
  const half = Math.atan2(dy, dx) / 2;
  return {
    data: [{x1: ndt + r / pxPerSec, y1: start,
            x2: ndt + (r * dx / len) / pxPerSec,
            y2: start + (r * dy / len) / pxPerEv}],
    bend: -half * 180 / Math.PI,
    // Out along the middle of the wedge, just clear of the arc.
    label: {x: ndt + (r + 14) * Math.cos(half) / pxPerSec,
            y: start + (r + 14) * Math.sin(half) / pxPerEv}
  };
}
```

``` js
// Mean decision time, over both responses, which is where the dashed line
// stands. With one boundary it is `boundary / drift`, the inverse Gaussian's
// mean, and the drift arrow lands exactly on it. With two it is shorter: the
// accumulations that would have taken longest are the ones that wander off and
// are absorbed at the other boundary, so they never contribute their long
// times to this mean. That is why the dashed line then stands to the left of
// the arrow - a mean path is no longer a mean time once the evidence can leave
// by the other door. `z` is the start point measured from the lower boundary,
// so this comes to `(edge / drift) * tanh(drift * edge)` at an unbiased start,
// and is shorter still at a biased one, whichever way it leans.
meanDT = {
  if (nbounds === 1) return boundary / drift;
  const z = boundary * bias;           // the start point, from the lower one
  return (boundary / drift) * (1 - Math.exp(-2 * drift * z)) /
                              (1 - Math.exp(-2 * drift * boundary)) - z / drift;
}
```

``` js
// The densities the paths are the first-passage times of, each drawn in
// evidence units on top of the boundary it belongs to.
//
// With one boundary that is the shifted Wald, i.e. dcogmod_invgaussian(t,
// drift, boundary, ndt). With two it is the Wiener first-passage density,
// dcogmod_ddm(t, drift, 2 * boundary, bias, ndt), defective at each boundary:
// each curve integrates to the share of the responses that boundary takes
// rather than to one. At an unbiased start the two have the same shape and
// differ only in area, in the ratio exp(drift * boundary), so everything the
// drift does to the choice is in the areas. A start point off the middle
// breaks that: the near boundary takes both more of the responses and the
// faster ones.
//
// The pair is pinned by its common maximum, so the taller curve always reaches
// `dheight` and the shorter one keeps its share. The vertical scale of a
// density means nothing on an evidence axis, and at drift 6, boundary 0.5 the
// peak is twelve times what it is at drift 1 - unpinned, it would leave the
// frame.
density = {
  const n = 500;
  const raw = [];
  let fmax = 0;
  for (let i = 0; i <= n; i++) {
    const t = cfg.tmax * i / n;
    const d = t - ndt;
    let up = 0, low = 0;
    if (d > 0) {
      if (nbounds === 1) {
        up = boundary / Math.sqrt(2 * Math.PI * Math.pow(d, 3)) *
             Math.exp(-Math.pow(boundary - drift * d, 2) / (2 * d));
      } else {
        // `boundary` is the separation here, which is the scale the whole
        // density is written in: what the two share is the standardised
        // density rescaled to it, with the drift's own factor. The upper
        // boundary is the lower boundary of the reflected process, which is
        // what flips the drift and the start point.
        const k = Math.exp(-drift * drift * d / 2) / (boundary * boundary);
        const u = d / (boundary * boundary);
        low = k * Math.exp(-drift * boundary * bias) * fpt0(u, bias);
        up = k * Math.exp(drift * boundary * (1 - bias)) * fpt0(u, 1 - bias);
      }
    }
    fmax = Math.max(fmax, up, low);
    raw.push({t: t, up: up, low: low});
  }
  const scale = fmax > 0 ? cfg.dheight / fmax : 0;
  return {
    up: raw.map(r => ({t: r.t, base: edge, y: edge + r.up * scale})),
    low: nbounds === 1 ? []
      : raw.map(r => ({t: r.t, base: -edge, y: -edge - r.low * scale}))
  };
}
```

``` js
// First-passage density at the lower boundary of the standardised Wiener
// process - no drift, unit separation, started at `w` - at time `u`. Every
// other diffusion's density is this one rescaled, which is the rescaling
// `density` above does. It is Navarro and Fuss (2009), the same pair of series
// R's `.ddm_lfpt0()` sums, but with fixed term counts rather than their
// error-driven ones: the figure never leaves 0 < u < 4, and over that stretch
// the counts below are well past the point where double precision stops
// noticing. The two series converge from opposite ends and both hold
// comfortably at u = 0.5, which is where the cheaper of them takes over.
fpt0 = function(u, w) {
  if (!(u > 0)) return 0;
  if (u < 0.5) {
    let s = 0;
    for (let k = -3; k <= 3; k++) {
      const wk = w + 2 * k;
      s += wk * Math.exp(-wk * wk / (2 * u));
    }
    return Math.max(0, s) / Math.sqrt(2 * Math.PI * u * u * u);
  }
  let s = 0;
  for (let k = 1; k <= 12; k++) {
    s += k * Math.exp(-k * k * Math.PI * Math.PI * u / 2) *
         Math.sin(k * Math.PI * w);
  }
  return Math.max(0, Math.PI * s);
}
```

``` js
// The Brownian increments the paths are built from, drawn once per volley and
// held while the sliders move. Two things follow from holding them. The same
// thirty trials are re-integrated under the new parameters rather than redrawn
// from fresh noise, so dragging `drift` shows this sample getting faster
// instead of a new sample that happens to be faster - which is the point of
// standing paths next to a density. And the expensive half of the work leaves
// the drag: a pointermove integrates a matrix that is already there.
//
// A row is as long as the longest run the sliders allow. The count below is
// `(tmax - ndt) / dt` and the ndt slider starts at 0, so the longest run is
// the whole of tmax and this is that same expression with ndt at its floor -
// exactly long enough, and only because the floor is 0. Move it and this has
// to move with it, or a path reads past the end of its row and integrates a
// NaN. The sqrt(dt) is in here rather than in the loop: these are increments
// of the process, not standard normals.
noise = {
  volley;                                       // a fresh draw with each volley
  const rnorm = d3.randomNormal(0, 1);
  const sd = Math.sqrt(cfg.dt);
  const nsteps = Math.ceil(cfg.tmax / cfg.dt);
  const step = () => sd * rnorm();
  return Array.from({length: cfg.npaths},
                    () => Float64Array.from({length: nsteps}, step));
}
```

``` js
// Euler-Maruyama paths for the same process, re-integrated whenever a slider
// moves. They set off from the start point, and a path ends at whichever
// boundary it touches first, carrying the height it ended at so that the dot
// which closes it knows which one to sit on.
trials = {
  const nsteps = Math.ceil((cfg.tmax - ndt) / cfg.dt);
  const floor = nbounds === 1 ? -Infinity : -edge;
  const rows = [];
  const hits = [];
  for (let k = 0; k < cfg.npaths; k++) {
    const dw = noise[k];
    let x = start;
    rows.push({id: k, t: 0, x: start});
    rows.push({id: k, t: ndt, x: start});
    for (let i = 1; i <= nsteps; i++) {
      x += drift * cfg.dt + dw[i - 1];
      const t = ndt + i * cfg.dt;
      if (x >= edge || x <= floor) {
        const end = x >= edge ? edge : floor;
        rows.push({id: k, t: t, x: end});
        hits.push({id: k, t: t, x: end});
        break;
      }
      rows.push({id: k, t: t, x: x});
    }
  }
  return {rows: rows, hits: hits};
}
```

## The Data

``` r

library(cogmod)
library(easystats)
library(ggplot2)
library(dplyr)
library(brms)
library(cmdstanr)

options(mc.cores = parallel::detectCores() - 2)
```

Decision making models jointly account for the **choice** that was made
and the **response time (RT)** it took to make it. Rather than
simulating data, we re-use the Wagenmakers et al. (2008) lexical
decision data (see the [RT-only
Models](https://dominiquemakowski.github.io/cogmod/articles/rt_models.md)
vignette) - but this time, instead of discarding the errors, we model
**choice** as *correct* vs. *error* responses. This is a common strategy
in the decision-making literature when the task itself does not have a
natural “left vs. right” stimulus category to map onto the two
accumulators/boundaries of the models below.

Evidence accumulation models are considerably more expensive to sample
than the RT-only models. We therefore use a smaller subset of the data
here, only three participant, so that the models below can be fit in a
reasonable amount of time for demonstration purposes.

``` r

set.seed(123)  # For reproducibility

# Experiment 1 of Wagenmakers et al. (2008), from rtdists. 
data(speed_acc, package = "rtdists")

df <- data.frame(
  Participant = as.integer(as.character(speed_acc$id)),
  Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[
    as.character(speed_acc$condition)]),
  RT = speed_acc$rt,
  Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat)),
  Frequency = unname(c(high = "High", low = "Low", very_low = "Very Low")[
    sub("^nw_", "", as.character(speed_acc$frequency))])
) 


df <- df[df$Participant %in% c(1, 2, 3) & df$RT <= 2, ]

# Show 10 first rows
head(df, 10)
#>    Participant Condition    RT Error Frequency
#> 1            1     Speed 0.700     0       Low
#> 2            1     Speed 0.392     1  Very Low
#> 3            1     Speed 0.460     0  Very Low
#> 4            1     Speed 0.455     0  Very Low
#> 5            1     Speed 0.505     1       Low
#> 6            1     Speed 0.773     0      High
#> 7            1     Speed 0.390     0      High
#> 8            1     Speed 0.587     1       Low
#> 9            1     Speed 0.603     0       Low
#> 10           1     Speed 0.435     0      High
```

``` r

ggplot(df, aes(x = RT, fill = Condition)) +
  geom_histogram(data = df[df$Error == 0, ], aes(y = after_stat(count) / (nrow(df) * 0.02)),
                 binwidth = 0.02, alpha = 0.6, position = "identity") +
  geom_histogram(data = df[df$Error == 1, ], aes(y = -after_stat(count) / (nrow(df) * 0.02)),
                 binwidth = 0.02, alpha = 0.6, position = "identity") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  labs(x = "RT (s)", y = "Distribution (Error - Correct)", fill = "Response") +
  scale_fill_manual(values = c("Accuracy"="#3F51B5", "Speed"="#F4511E")) +
  theme_minimal()
```

![](decision_making_files/figure-html/unnamed-chunk-30-1.png)

Errors are much rarer than correct responses (especially in the
`Accuracy` condition), which can be problematic for accurate
estimations.

## Models

All five models below use `dec(Error)` to indicate the two-choice
outcome (`0` = Correct, `1` = Error), and all five share the
parameterization used by the RT-only families (see
`vignette("outliers")`): `ndt` is estimated directly, in seconds, and a
`poutlier` parameter mixes in an outlier process so that `ndt` is not
capped at the fastest observed response. The outlier component is a half
Normal with a fixed scale of 0.2 s, so **reaction times must be in
seconds**.

[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
and
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
read the family off the formula, so the same three lines set up any of
them. Both are worth using rather than hand-written priors and
`init = 0`: on a log link `init = 0` starts `ndt` at `exp(0) = 1`
second, above nearly every observed RT, and each of these families has
at least one direction the likelihood is flat in that `brms` would
otherwise leave improper.

### Drift Diffusion Model (DDM)

The DDM assumes that evidence accumulates towards one of two boundaries
at a rate `mu` (drift rate). `boundary` is the boundary separation
(higher = more cautious), `bias` is the starting point between the two
boundaries (`0.5` = unbiased), and `ndt` is the non-decision time, in
seconds.

Two conventions invert the signs one might expect. Following `brms`’s
own
[`wiener()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family, the response coded `1` in `dec()` - here the **error** - is the
*upper* boundary, and `bias` is measured from the *lower* one. Good
performance therefore shows up as a **negative** `mu`, and `bias > 0.5`
puts the start point closer to the error boundary, making errors
*faster* than correct responses (`bias < 0.5` makes them slower; at
exactly `0.5` the two conditional RT distributions are identical).

Note that we use the “simple” 4-parameter DDM here, which does not
include between-trial variability in drift rate, starting point, or
non-decision time, hence the `sigmadrift`, `sigmabias`, and `sigmandt`
parameters are fixed to `0`. Estimating these parameters is possible,
but considerably more expensive and often unnecessary for many
applications. Note that `sigmandt` is the between-trial *range* of the
non-decision time in seconds (`st0`), with `ndt` its lower bound.

``` r

f <- bf(
  RT | dec(Error) ~ Condition,
  boundary ~ Condition,
  bias ~ 1,
  ndt ~ 1,
  sigmadrift = 0,
  sigmabias = 0,
  sigmandt = 0,
  family = cogmod_ddm()
)

m_ddm <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_ddm <- brms::add_criterion(m_ddm, "loo")  # Add model performance criterion
```

### DDM with Drift Variability (DDM-5)

The three variability parameters each produces a specific effect on the
*relative* speed of correct and error responses. Between-trial
variability in the **drift rate** (`sigmadrift`, Ratcliff’s $`s_v`$)
makes errors *slower* than correct responses, because errors are then
contributed disproportionately by the trials that happened to draw a low
drift. Variability in the **starting point** (`sigmabias`, $`s_z`$) does
the opposite, producing *faster* errors, and variability in the
**non-decision time** (`sigmandt`, $`s_{t0}`$) mostly affects the
leading edge of the distribution. Which one to free is therefore an
empirical question with a visible answer, and here the errors are
slightly *slower* than the correct responses - so `sigmadrift` is the
parameter to free.

We keep `sigmabias` and `sigmandt` fixed at `0`, for two reasons.
Statistically they are weakly identified with only few error trials, and
computationally the exact zero matters: `cogmod`’s Stan code falls back
to the dedicated (and much cheaper) drift-variability-only density when
both are exactly `0`, and to adaptive numerical quadrature otherwise. A
prior *concentrated near* zero would pay the full cost of the
7-parameter form without buying anything.

``` r

f <- bf(
  RT | dec(Error) ~ Condition,
  boundary ~ Condition,
  bias ~ 1,
  ndt ~ 1,
  sigmadrift ~ 1,
  sigmabias = 0,
  sigmandt = 0,
  family = cogmod_ddm()
)

# cogmod_priors() already supplies normal(0, 1) for sigmadrift, for the same
# reason it fences the other flat directions.
m_ddm5 <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_ddm5 <- brms::add_criterion(m_ddm5, "loo")  # Add model performance criterion
```

### LogNormal Race (LNR)

The LNR is a somewhat simpler model. It is similar to the LBA, but each
accumulator’s finishing time is drawn directly from a LogNormal
distribution instead of a ballistic accumulation process. `mu`
(`nuzero`) and `nuone` are the (inverse log-space mean) processing
speeds for the “Error” and “Correct” accumulators, and
`sigmazero`/`sigmaone` their log-space SDs.

The LNR also has a `sigmabias` parameter, a Uniform start-point range
below a threshold pinned one unit above it, which turns it into the LBA
with LogNormal drift rates ([Heathcote & Love,
2012](https://doi.org/10.3389/fpsyg.2012.00292)). It is fixed at zero
here, which is the LNR proper: a start-point range is hard to identify
from the shape of the RT distribution alone and adds a flat direction to
the likelihood at zero, so it is worth estimating only with a lot of
data or a specific hypothesis about start-point variability (see
[`?cogmod_lnr`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)).

Like every family in this package, the LNR is fit with `ndt` and
`poutlier`: `ndt` is estimated in seconds, with no upper bound tied to
the fastest observed response, and `poutlier` is the proportion of
trials attributed to a contaminant guessing process (see
`vignette("outliers")`).
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
and
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
both read the family off `f`, so they need no family-specific setup -
and `init = 0` is actively harmful here: on the log link it starts `ndt`
at `exp(0) = 1` second, above nearly every observed RT.

``` r

f <- bf(
  RT | dec(Error) ~ Condition,
  nuone ~ Condition,
  sigmazero ~ 1,
  sigmaone ~ 1,
  sigmabias = 0,
  ndt ~ Condition,
  family = cogmod_lnr()
)

m_lnr <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_lnr <- brms::add_criterion(m_lnr, "loo")  # Add model performance criterion
```

### Linear Ballistic Accumulator (LBA)

The LBA assumes two independent accumulators (one per choice) that race
towards a common threshold `b` (`sigmabias` = start-point range `A`,
`boundary` = extra distance so that `b = A + boundary`). `mu` and
`driftone` are the mean drift rates for the “Correct” and “Error”
accumulators, and `sigmazero`/ `sigmaone` their between-trial drift
variability.

Note that `sigmazero` is **fixed to 1** below. The evidence scale of an
LBA is arbitrary - multiply the drifts, their SDs, the start-point range
and the threshold by any constant and every finishing time is
unchanged - so the six parameters are identified only up to a common
factor. Priors make the posterior proper; only fixing one of them
identifies the scale.

Two more things are worth knowing before reading the estimates. First,
each drift rate is a Normal **truncated at zero** (the convention of
`rtdists`, `DMC` and `EMC2`; see
[`?rcogmod_lba2`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)),
and for an accumulator that rarely wins - the error accumulator here,
with a 5% error rate - the truncated Normal’s location and scale are
identified only through their ratio. `driftone` will therefore come out
well below zero with a wide `sigmaone`, and the pair should be read
together, as the shape of a distribution of small positive rates, rather
than `driftone` alone as a mean drift. Fixing `sigmaone = 1` as well is
sometimes worth doing.

Second, in our case, the `boundary` and `ndt` parameter are strongly
correlated, which makes sampling difficult and slow. In these cases,
setting `metric = "dense_e"` can help (x2 speed-ups in our case). This
setting is worth trying in other models as well. While it might pay off
for low-dimensional posteriors with strong correlations, for a hierarchy
with hundreds of participant-level parameters the dense matrix has more
entries to estimate than warmup can pin down, and the default might be
the safer choice. The
[performance](https://dominiquemakowski.github.io/cogmod/articles/performance.md)
article explains the trade-off and reports what it bought across the
families in a local benchmark.

``` r

f <- bf(
  RT | dec(Error) ~ Condition,
  driftone ~ Condition,
  sigmazero = 1,
  sigmaone ~ 1,
  sigmabias ~ 1,
  boundary ~ 1,
  ndt ~ 1,
  family = cogmod_lba2()
)

m_lba <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 500, backend = "cmdstanr",
  metric = "dense_e"  # see above
)

m_lba <- brms::add_criterion(m_lba, "loo")  # Add model performance criterion
```

### Racing Diffusion Model (RDM)

The RDM is the LBA’s stochastic counterpart - each accumulator
integrates evidence through the DDM’s random walk process instead of
accumulating linearly. It otherwise keeps the racing architecture: two
independent accumulators, a common threshold `b`, and a start point
drawn from `Uniform(0, sigmabias)`. The consequence of swapping the
ballistic path for a diffusing one is that the noise now lives *within*
a trial rather than between trials. That is the point of the model:
Tillman et al. ([2020](https://doi.org/10.3758/s13423-020-01738-8)) show
that within-trial variability alone accounts for the benchmark choice-RT
phenomena, without the between-trial drift variability that the LBA
(`sigmazero`/`sigmaone`) and the full DDM (`sigmadrift`) need. It is
therefore the more parsimonious race: `mu` and `driftone` are the drift
rates for the “Correct” and “Error” accumulators, and there is no drift
variability parameter to estimate.

Because each accumulator is a Wald (shifted inverse Gaussian) process,
the drift rates use a softplus link and are constrained to be
non-negative. Like the LNR, `ndt` is estimated in seconds, and
`poutlier` is the proportion of trials attributed to a contaminant early
responses (see `vignette("outliers")`).

``` r

f <- bf(
  RT | dec(Error) ~ Condition,
  driftone ~ Condition,
  sigmabias ~ 1,
  boundary ~ 1,
  ndt ~ 1,
  family = cogmod_rdm()
)

m_rdm <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_rdm <- brms::add_criterion(m_rdm, "loo")  # Add model performance criterion
```

**Identifiability of `sigmabias` and `boundary`**: This concerns the LBA
as much as the RDM, since both are parameterized the same way. The two
parameters enter the model only through the threshold
`b = boundary + sigmabias`, so they trade off almost freely, potentially
ruining model convergence.

This is a different identifiability problem from the model’s overall
*scale* invariance - multiplying every drift, its SD, the start-point
range and the threshold by a constant leaves every finishing time
unchanged - which the LBA literature resolves by fixing one drift-rate
SD to `1` (Brown & Heathcote,
[2008](https://doi.org/10.1016/j.cogpsych.2007.12.002)), the convention
already used for `sigmazero` above. There is no equivalent standard fix
for the `boundary`/`sigmabias` split itself. Toolboxes that estimate by
maximum likelihood (e.g. `rtdists`) just let both float freely and
accept the estimation noise that comes with it; Bayesian hierarchical
packages built around this model class (e.g. `EMC2`, the successor to
`DMC`) do not fix either parameter outright either - they regularize the
ridge with weakly-informative priors and partial pooling across
participants and conditions rather than pin the split down directly.
`cogmod` takes the same route at the level of a single fit: a prior on
`sigmabias` might help, and for the RDM
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
supplies a weakly informative one (`normal(0, 1)`) on the softplus
scale.

A prior of that shape does not make `sigmabias` itself trustworthy,
though, so prefer `boundary + sigmabias` whenever you interpret a
threshold or compare one across conditions, and treat the split between
the two as weakly-determined.

**When errors are too few to identify the second accumulator**: every
race above has to estimate an error process, and when errors are scarce
that process is informed by almost nothing. On the Accuracy condition of
these data, which carries about 75 errors, an RDM left to itself runs
`driftone` down to the floor of its link (a drift of 0.0001) with half
of the transitions divergent;
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
fences that with a prior on the error drift (see
[`?rcogmod_rdm`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)).
The structural alternative is to stop modelling the error process
altogether and fit the RT-only family with the errors as
**right-censored** correct responses, `bf(RT | cens(Error) ~ ...)`: an
error then says only that the correct process had not finished yet, so
there is no error accumulator to run away. On
[`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
this is the censored shifted Wald of Miller et
al. ([2018](https://doi.org/10.1177/0146621617710465)). It is described
in the *Censored Shifted Wald* section of the [RT-only
Models](https://dominiquemakowski.github.io/cogmod/articles/rt_models.md)
vignette, together with the one check to run first: censoring can only
produce errors *slower* than correct responses, so where errors are
faster - a low boundary, a biased start point - stay with the race.

## Model Comparison

> Everything below is fitted to a subset of data, with no random
> effects, short chains, and predictors on only some parameters - with
> the goal of demonstrating the workflow. The “findings” are thus not to
> be taken at face value.

### Model Fit

``` r

loo::loo_compare(m_ddm, m_ddm5, m_lba, m_lnr, m_rdm) |>
  parameters(include_ENP = TRUE)
#> # Fixed Effects
#> 
#> Name   |   LOOIC |   ENP |    ELPD | Difference | Difference_SE |      p
#> ------------------------------------------------------------------------
#> m_lba  | -2657.6 |  7.38 | 1328.82 |       0.00 |          0.00 |       
#> m_lnr  | -2632.2 | 10.41 | 1316.08 |     -12.74 |         10.90 | 0.242 
#> m_ddm5 | -2506.6 | 11.06 | 1253.31 |     -75.52 |         18.84 | < .001
#> m_ddm  | -2436.2 | 10.66 | 1218.08 |    -110.75 |         21.43 | < .001
#> m_rdm  | -2419.6 |  6.92 | 1209.82 |    -119.00 |         18.14 | < .001
```

### Sampling Duration

Choice+RT models are considerably more expensive to sample than RT-only
models: the DDM relies on Stan’s `wiener_lpdf`, which is comparatively
slow, while the LBA, LNR and RDM likelihoods involve evaluating both a
“winner” density and a “loser” survival function for every observation.
Among the three races, the LNR is by far the cheapest, since its density
and survival are just LogNormal ones, while the RDM and the LBA are much
slower. Freeing `sigmadrift` is not free either: DDM-5 costs is
significantly slower than the simple DDM, though it remains cheaper than
either of the two slow races.

``` r

models <- list(
  DDM = m_ddm, `DDM-5` = m_ddm5, LNR = m_lnr, LBA = m_lba, RDM = m_rdm
)
model_levels <- names(models)

duration <- do.call(rbind, lapply(model_levels, function(nm) {
  data_modify(attributes(models[[nm]]$fit)$metadata$time$chain, Model = nm)
})) |>
  data_modify(Model = factor(Model, levels = model_levels), Minutes = total / 60)

duration_range <- duration |>
  summarize(
    duration_min = min(Minutes),
    duration_median = median(Minutes),
    duration_max = max(Minutes),
    .by = Model
  )

quality <- do.call(rbind, lapply(model_levels, function(nm) {
  est <- models[[nm]]$criteria$loo$estimates
  data.frame(Model = nm, elpd = est["elpd_loo", "Estimate"], elpd_se = est["elpd_loo", "SE"])
})) |>
  data_modify(Model = factor(Model, levels = model_levels))

fit_summary <- merge(duration_range, quality, by = "Model")

fit_summary |>
  ggplot(aes(x = duration_median, y = elpd, color = Model)) +
  geom_errorbar(aes(xmin = duration_min, xmax = duration_max), orientation = "y") +
  geom_errorbar(aes(ymin = elpd - elpd_se, ymax = elpd + elpd_se), width = 0) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = Model), size = 3.2, show.legend = FALSE) +
  scale_color_material_d(guide = "none") +
  labs(
    x = "Sampling Duration per Chain (min) - median, range across the 4 chains",
    y = "Fit Quality (elpd_loo ± 1 SE)"
  ) +
  theme_minimal()
```

![](decision_making_files/figure-html/unnamed-chunk-43-1.png)

### Posterior Predictive Check

Each model predicts a **pair** of outcomes per trial, so
`estimate_prediction()` returns them in a `Component` column (`"rt"` and
`"response"`), which we pivot back into two columns.

Correct responses are drawn upwards and errors downwards, with the
observed data as histograms and the five models as overlaid density
lines (faint: one per posterior draw; bold: pooled over draws). Rather
than normalizing each half separately, both are scaled by the
**proportion** of that response, so that the area under each curve
equals its predicted frequency - which is why the error half is roughly
a twelfth of the size of the correct one. This makes the plot a check on
the *joint* distribution of choices and RTs: a model can only match it
by getting the error rate *and* the shape of both RT distributions
right.

Code

``` r

pred <- rbind(
  estimate_prediction(m_ddm, data = df, iterations = 50, 
                      keep_iterations = TRUE, ci = NULL) |>
    reshape_iterations() |>
    data_modify(Model = "DDM"),
  estimate_prediction(m_ddm5, data = df, iterations = 50, 
                      keep_iterations = TRUE, ci = NULL) |>
    reshape_iterations() |>
    data_modify(Model = "DDM-5"),
  estimate_prediction(m_lnr, data = df, iterations = 50, 
                      keep_iterations = TRUE, ci = NULL) |>
    reshape_iterations() |>
    data_modify(Model = "LNR"),
  estimate_prediction(m_lba, data = df, iterations = 50, 
                      keep_iterations = TRUE, ci = NULL) |>
    reshape_iterations() |>
    data_modify(Model = "LBA"),
  estimate_prediction(m_rdm, data = df, iterations = 50, 
                      keep_iterations = TRUE, ci = NULL) |>
    reshape_iterations() |>
    data_modify(Model = "RDM")
) |>
  datawizard::data_select(select = c("Row", "Component", "Condition", "iter_value", "iter_group", "iter_index", "Model")) |>
  datawizard::data_to_wide(id_cols = c("Row", "iter_group", "Model", "Condition"), values_from = "iter_value", names_from = "Component") |> 
  data_modify(Model = factor(Model, levels = c("DDM", "DDM-5", "LNR", "LBA", "RDM")))

# The LBA occasionally predicts enormous RTs (a near-zero drift rate takes a
# very long time to reach the threshold). `stat_density()` spreads its
# evaluation grid over the whole x-axis range, so a single extreme draw would
# flatten every curve.
pred <- data_filter(pred, rt < 3)

correct <- pred[pred$response == 0, ]
error <- pred[pred$response == 1, ]

n_obs <- nrow(df)  # Trials per posterior draw
n_iter <- length(unique(pred$iter_group))
bw <- 0.02  # Histogram bin width

# Dividing the counts by the *total* number of trials (rather than by the
# number of trials of that response) scales each half by its own frequency.
p <- ggplot(df, aes(x = RT)) +
  # Observed data
  geom_histogram(data = df[df$Error == 0, ], aes(y = after_stat(count) / (nrow(df) * 0.02),
                                                 fill = Condition),
                 binwidth = 0.02, alpha = 0.6, position = "identity") +
  geom_histogram(data = df[df$Error == 1, ], aes(y = -after_stat(count) / (nrow(df) * 0.02),
                                                 fill = Condition),
                 binwidth = 0.02, alpha = 0.6, position = "identity") +
  # One faint line per posterior draw
  # geom_line(data = correct,
  #           aes(x = rt, y = after_stat(count) / n_obs, color = Model,
  #               group = interaction(Model, iter_group, Condition)),
  #           stat = "density", alpha = 0.05, linewidth = 0.3) +
  # geom_line(data = error,
  #           aes(x = rt, y = -after_stat(count) / n_obs, color = Model,
  #               group = interaction(Model, iter_group, Condition)),
  #           stat = "density", alpha = 0.05, linewidth = 0.3) +
  # Posterior predictive density, pooled over draws
  geom_line(data = correct,
            aes(x = rt, y = after_stat(count) / (n_obs * n_iter), 
                color = Model, linetype = Condition),
            stat = "density", linewidth = 1.2) +
  geom_line(data = error,
            aes(x = rt, y = -after_stat(count) / (n_obs * n_iter), 
                color = Model, linetype = Condition),
            stat = "density", linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_color_material_d(palette = "rainbow") +
  scale_fill_manual(values = c("Accuracy"="#3F51B5", "Speed"="#F4511E")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2)),
         linetype = "none") +
  coord_cartesian(xlim = c(0.25, 1.25)) +
  labs(x = "RT (s)", y = "Density (up = Correct, down = Error)", color = "Model") +
  # facet_wrap(~Model) +
  theme_minimal()
p
```

![](../reference/figures/decision_making1.png)

The mismatch discussed above is visible in the lower half: the **DDM**
(purple) places its error density well to the left of the observed
errors, while the three races put theirs on top of them. **DDM-5**
(blue) sits between the two - it recovers the average error timing but
spreads the errors too widely. In the upper half all five are close,
which is the point worth taking away: a model’s problem can be invisible
if one only looks at the RTs of correct responses.
