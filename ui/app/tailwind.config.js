/** @type {import('tailwindcss').Config} */
module.exports = {
  prefix: 'tw-',
  purge: ['./src/**/*.vue'],
  content: ['./index.html', './src/**/*.{vue,js,ts,jsx,tsx}'],
  darkMode: false, // or 'media' or 'class'
  theme: {
    extend: {
      maxWidth: {
        'screen-xxl': '1600px',
      },
    },
  },
  variants: {},
  plugins: [],
}
