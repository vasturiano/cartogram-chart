import { select as d3Select, event as d3Event } from 'd3-selection';
import 'd3-transition';
import { geoMercator } from 'd3-geo';
import cartogramModule from './lib/topogram/d3-cartogram.js';
import Kapsule from 'kapsule';
import accessorFn from 'accessor-fn';

const d3Cartogram = cartogramModule.cartogram; // Unwrap CJS from ES module

const ANIMATION_DURATION = 1200;

export default new Kapsule({

  props: {
    width: { default: window.innerWidth },
    height: { default: window.innerHeight },
    iterations: { default: 20 },
    projection: { default: geoMercator()
      .scale(Math.min((window.innerWidth - 3) / (2 * Math.PI), (window.innerHeight - 3) / (1.2 * Math.PI)))
      .translate([window.innerWidth / 2, window.innerHeight / 1.5])
    },
    topoJson: {},
    topoObjectName: {},
    value: { default: 1 },
    color: { default: 'lightgrey' },
    label: { default: '' },
    valFormatter: { default: n => n },
    units: { default: '' },
    tooltipContent: { default: d => '' },
    onClick: { default: d => {} }
  },

  init(domNode, state) {
    state.cartogram = d3Cartogram()
      .properties(d => d.properties);

    // Dom
    state.svg = d3Select(domNode).append('svg')
      .attr('class', 'cartogram');

    // tooltips
    state.tooltip = d3Select('body')
      .append('div')
      .attr('class', 'chart-tooltip cartogram-tooltip');

    // tooltip cleanup on unmount
    domNode.addEventListener ('DOMNodeRemoved', function(e) {
      if (e.target === this) { state.tooltip.remove(); }
    });

    state.svg.on('mousemove', () => {
      state.tooltip
        .style('left', d3Event.pageX + 'px')
        .style('top', d3Event.pageY + 'px');
    });
  },

  update(state) {
    const valueOf = accessorFn(state.value);
    const colorOf = accessorFn(state.color);
    const labelOf = accessorFn(state.label);
    const tooltipContentOf = accessorFn(state.tooltipContent);

    state.svg
      .attr('width', state.width)
      .attr('height', state.height);

    if (!state.topoJson) return; // No features to render

    const topoObject = state.topoJson.objects[state.topoObjectName] || Object.values(state.topoJson.objects)[0];
    if (!topoObject) {
      console.warn('Unable to find topology object in TopoJson');
      return;
    }

    state.cartogram
      .projection(state.projection)
      .value(valueOf);

    const features = state.svg.selectAll('path.feature')
      .data(state.cartogram
        .iterations(1) // Initialize new features non-distorted
        (state.topoJson, topoObject.geometries).features
      );

    features.exit().remove();

    const newFeatures = features.enter().append('path')
      .attr('class', 'feature')
      .style('fill', 'lightgrey')
      .attr('d', state.cartogram.path)
      .on('mouseover', (feature) => {
        const label = labelOf(feature);
        const extraContent = tooltipContentOf(feature);
        state.tooltip.style('display', 'inline');
        state.tooltip.html(`
          ${label ? `<b>${label}</b>:` : ''}
          ${state.valFormatter(valueOf(feature))}
          ${state.units}
          ${extraContent ? `<br/><br/>${extraContent}` : ''}
        `);
      })
      .on('mouseout', () => { state.tooltip.style('display', 'none'); })
      .on('click', state.onClick);

    features.merge(newFeatures)
      .data(state.cartogram
        .iterations(state.iterations) // distort all features
        (state.topoJson, topoObject.geometries).features
      )
      .transition().duration(ANIMATION_DURATION)
        .style('fill', colorOf)
        .attr('d', state.cartogram.path);
  }

});