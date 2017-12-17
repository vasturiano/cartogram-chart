# Cartogram Chart

<p align="center">
     <a href="https://vasturiano.github.io/cartogram-chart/example/world-population/"><img width="70%" src="https://vasturiano.github.io/cartogram-chart/example/world-population/screenshot.png"></a>
</p>

An interactive contiguous cartogram reusable chart for visualizing geographical data. 
Given a [TopoJson](https://github.com/topojson/topojson/wiki) topology, renders its shapes with distorted areas according to a value associated with each feature. The geo projection to be used is configurable using the `projection` property.
Uses [topogram](https://github.com/shawnbot/topogram) for the layout computation.

[![NPM](https://nodei.co/npm/cartogram-chart.png?compact=true)](https://nodei.co/npm/cartogram-chart/)

## Quick start

```
import Cartogram from 'cartogram-chart';
```
or
```
Cartogram = require('cartogram-chart');
```
or even
```
<script src="//unpkg.com/cartogram-chart"></script>
```
then
```
const myChart = Cartogram();
myChart
    .data(<myData>)
    (<myDOMElement>);
```

## API reference

| Method | Description | Default |
| --- | --- | --- |
| <b>width</b>([<i>number</i>]) | Getter/setter for the chart width in px. | *&lt;window width&gt;* |
| <b>height</b>([<i>number</i>]) | Getter/setter for the chart height in px. | *&lt;window height&gt;* |
| <b>topoJson</b>([<i>object</i>]) | Getter/setter for the [TopoJson](https://github.com/topojson/topojson/wiki) topology. Without this property no shapes are rendered. | |
| <b>topoObjectName</b>([<i>string</i>]) | Getter/setter for the object name in the `topoJson.objects` structure to use. | *&lt;first object&gt;* |
| <b>projection</b>([<i>object</i>]) | Getter/setter for the geographical [projection](https://github.com/d3/d3-geo-projection) to use for rendering. | `geoMercator` <i>(centered on prime meridian, slightly tilted towards the northern hemisphere)</i> |
| <b>iterations</b>([<i>number</i>]) | Getter/setter for the number of iterations to run the algorithm for. Higher numbers distorts the areas closer to their associated value, at the cost of performance. | 20 |
| <b>value</b>([<i>number</i>, <i>string</i> or <i>fn</i>]) | Getter/setter for a feature's value accessor. The shape area size is distorted according to this property. Supports either a fixed <i>numeric value</i>, a <i>string</i> indicating the features's object attribute, or a `function(feature)` which should return a numeric value. | 1 |
| <b>color</b>([<i>string</i> or <i>fn</i>]) | Getter/setter for a feature's color accessor, used to color the shapes. | `lightgrey` |
| <b>label</b>([<i>string</i> or <i>fn</i>]) | Getter/setter for a feature's label accessor, used to display a shape's name on its tooltip. | |
| <b>valFormatter</b>([<i>function</i>]) | Getter/setter for the number formatter `function(n)`, to show values in the tooltip. | `n => n` |
| <b>units</b>([<i>string</i>]) | Getter/setter for the value units, to include in the tooltip. | |
| <b>tooltipContent</b>([<i>string</i> or <i>fn</i>]) | Getter/setter for a feature's tooltip content accessor. Use this to specify extra content in each of the shape's tooltips in addition to the feature label and value that is included by default. | |
| <b>onClick</b>([<i>function</i>]) | Getter/setter for the callback `function(feature)` to trigger when clicking on a shape. | - |