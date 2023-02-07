import resolve from '@rollup/plugin-node-resolve';
import commonJs from '@rollup/plugin-commonjs';
import babel from '@rollup/plugin-babel';
import postCss from 'rollup-plugin-postcss';
import terser from '@rollup/plugin-terser';

import pkg from './package.json' assert { type: 'json' };
const { name, homepage, version, dependencies } = pkg;

const umdConf = {
  format: 'umd',
  name: 'Cartogram',
  banner: `// Version ${version} ${name} - ${homepage}`
};

export default [
  { // UMD
    input: 'src/index.js',
    output: [
      {
        ...umdConf,
        file: `dist/${name}.js`,
        sourcemap: true
      },
      { // minify
        ...umdConf,
        file: `dist/${name}.min.js`,
        plugins: [terser({
          output: { comments: '/Version/' }
        })]
      }
    ],
    plugins: [
      postCss({ plugins: [] }),
      babel({ exclude: 'node_modules/**' }),
      resolve(),
      commonJs()
    ]
  },
  { // ES module
    input: 'src/index.js',
    output: [
      {
        format: 'es',
        file: `dist/${name}.mjs`
      }
    ],
    external: [
      ...Object.keys(dependencies || {}),
      'three/examples/jsm/controls/DragControls.js'
    ],
    plugins: [
      postCss({ plugins: [] }),
      babel()
    ]
  }
];