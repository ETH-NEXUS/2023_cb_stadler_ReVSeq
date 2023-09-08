<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'

const props = defineProps({
  domain: {
    type: String,
    required: true,
  },
})

const {t} = useI18n()

const coreStore = useCoreStore()

const {selected_plate, selected_sample} = storeToRefs(coreStore)

const downloadFile = async (filePath: string) => {
  await coreStore.downloadFile(filePath)
}
</script>

<template>
  <ul v-if="props.domain === 'plate' && selected_plate">
    <li :key="file.path" v-for="file in selected_plate.files" class="q-my-sm">
      <button class="link-file" style="cursor: pointer" @click="downloadFile(file.path)">
        {{ file.path.split('/').pop() }}
      </button>
    </li>
  </ul>
  <ul v-if="props.domain === 'sample' && selected_sample">
    <li :key="file.path" v-for="file in selected_sample[0].files" class="q-my-sm">
      <button class="link-file" style="cursor: pointer" @click="downloadFile(file.path)">
        {{ file.path.split('/').pop() }}
      </button>
    </li>
  </ul>
</template>

<style>
.link-file {
  cursor: pointer !important;
  color: #328ee7;
}

.link-file:hover {
  cursor: pointer !important;
  color: #03396c;
}
</style>
