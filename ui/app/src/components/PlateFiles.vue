<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'
import {ref} from 'vue'
import {api} from 'src/boot/axios'

const props = defineProps({
  domain: {
    type: String,
    required: true,
  },
})

const {t} = useI18n()

const coreStore = useCoreStore()

const {selected_plate, selected_sample} = storeToRefs(coreStore)
const progress = ref(0);

// const downloadFile = async (filePath: string) => {
//   await coreStore.downloadFile(filePath)
// }


const downloadFile = async (path: string) => {
  try {
    progress.value = 0;
    const encodedPath = encodeURIComponent(path)
    const res = await api.get(`/api/download/${encodedPath}/`, {
      responseType: 'blob',
        onDownloadProgress: progressEvent => {
        if (progressEvent.total) {
          progress.value = progressEvent.loaded / progressEvent.total;
        } else {
          progress.value = 0.1;
        }
      }
    })

       const url = window.URL.createObjectURL(new Blob([res.data]))
      const link = document.createElement('a')
      link.href = url
      link.setAttribute('download', path.split('/').pop() || '')
      document.body.appendChild(link)
      link.click()
  } catch (error) {
    console.error('Error downloading the file:', error)
  } finally {
    progress.value = 0; // Reset progress after download completes or fails
  }
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
    <q-ajax-bar :value="progress"  position="bottom"
      color="accent"
      size="10px"
      />
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
