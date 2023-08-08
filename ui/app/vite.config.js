import { defineConfig } from 'vite';
import vue from '@vitejs/plugin-vue';
import { quasar } from '@quasar/vite-plugin';
import tailwindcss from 'tailwindcss';
import autoprefixer from 'autoprefixer';

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue(), quasar()],
  css: {
    postcss: {
      plugins: [tailwindcss('./tailwind.config.js'), autoprefixer],
    },
  },
});
