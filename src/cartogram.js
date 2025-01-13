import { select as d3Select } from 'd3-selection';
import 'd3-transition'; // extends d3-selection prototype
import { geoMercator } from 'd3-geo';
import { cartogram as d3Cartogram } from 'topogram';
import Kapsule from 'kapsule';
import accessorFn from 'accessor-fn';
import Tooltip from 'float-tooltip';

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
    label: { triggerUpdate: false },
    valFormatter: { default: n => n, triggerUpdate: false },
    units: { default: '', triggerUpdate: false },
    tooltipContent: { triggerUpdate: false },
    onClick: {}
  },

  init(domNode, state) {
    state.cartogram = d3Cartogram()
      .properties(d => d.properties);

    // Dom
    const el = d3Select(domNode)
      .append('div').attr('class', 'cartogram');

    state.svg = el.append('svg');

    // tooltip
    state.tooltip = new Tooltip(el);
  },

  update(state) {
    const valueOf = accessorFn(state.value);
    const colorOf = accessorFn(state.color);

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
      .on('mouseover', (ev, feature) => {
        const valueOf = accessorFn(state.value);
        const labelOf = accessorFn(state.label);
        const tooltipContentOf = accessorFn(state.tooltipContent);

        const label = labelOf(feature);
        const extraContent = tooltipContentOf(feature);
        state.tooltip.content(!label && !extraContent ? null : `
          ${label ? `<b>${label}</b>:` : ''}
          ${state.valFormatter(valueOf(feature))}
          ${state.units}
          ${extraContent ? `<br/><br/>${extraContent}` : ''}
        `);
      })
      .on('mouseout', () => { state.tooltip.content(null); })
      .on('click', (ev, d) => state.onClick && state.onClick(d));

    features.merge(newFeatures)
      .data(state.cartogram
        .iterations(state.iterations) // distort all features
        (state.topoJson, topoObject.geometries).features
      )
      .style('cursor', state.onClick ? 'pointer' : null)
      .transition().duration(ANIMATION_DURATION)
        .style('fill', colorOf)
        .attr('d', state.cartogram.path);
  }

});