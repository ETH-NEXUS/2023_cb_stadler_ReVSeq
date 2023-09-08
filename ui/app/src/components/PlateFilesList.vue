<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'

import {watchEffect} from 'vue'

const {t} = useI18n()

const coreStore = useCoreStore()

const {selected_plate} = storeToRefs(coreStore)
watchEffect(() => {
  console.log('Reactive property changed:', selected_plate.value)
})

const downloadFile = async () => {
  alert('download file')
}
</script>

<template>
  <div class="link-file" @click="downloadFile">jjjjj</div>
  <ul>
    <li :key="file.path" v-for="file in selected_plate.files" class="q-my-sm">
      <button class="link-file" style="cursor: pointer" @click="downloadFile">{{ file.path }}</button>
    </li>
  </ul>
</template>

<style>
@import 'tailwindcss/base.css';
@import 'tailwindcss/components.css';
@import 'tailwindcss/utilities.css';
.link-file {
  cursor: pointer !important;
  color: #328ee7;
}

.link-file:hover {
  cursor: pointer !important;
  color: #03396c;
}
</style>
