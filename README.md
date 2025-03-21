Cartogram Chart
==============

[![NPM package][npm-img]][npm-url]
[![Build Size][build-size-img]][build-size-url]
[![NPM Downloads][npm-downloads-img]][npm-downloads-url]

<p align="center">
     <a href="https://vasturiano.github.io/cartogram-chart/example/world-population/"><img width="70%" src="https://vasturiano.github.io/cartogram-chart/example/world-population/screenshot.png"></a>
</p>

An interactive contiguous cartogram reusable chart for visualizing geographical data. 

Given a [TopoJson](https://github.com/topojson/topojson/wiki) topology, renders its shapes with distorted areas according to a value associated with each feature. The geo projection to be used is configurable using the `projection` property.

Uses [Shawn Allen](https://github.com/shawnbot)'s [topogram](https://github.com/shawnbot/topogram) for the algorithm computation.

## Quick start

```js
import Cartogram from 'cartogram-chart';
```
or using a *script* tag
```html
<script src="//cdn.jsdelivr.net/npm/cartogram-chart"></script>
```
then
```js
const myChart = new Cartogram(<myDOMElement>)
  .topoJson(<myTopology>);
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

## Giving Back

[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donate_SM.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&business=L398E7PKP47E8&currency_code=USD&source=url) If this project has helped you and you'd like to contribute back, you can always [buy me a ☕](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&business=L398E7PKP47E8&currency_code=USD&source=url)!


[npm-img]: https://img.shields.io/npm/v/cartogram-chart
[npm-url]: https://npmjs.org/package/cartogram-chart
[build-size-img]: https://img.shields.io/bundlephobia/minzip/cartogram-chart
[build-size-url]: https://bundlephobia.com/result?p=cartogram-chart
[npm-downloads-img]: https://img.shields.io/npm/dt/cartogram-chart
[npm-downloads-url]: https://www.npmtrends.com/cartogram-chart
