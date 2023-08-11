<script setup lang="ts">
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'
import {ref, onMounted} from 'vue'

onMounted(async () => {
  await coreStore.getPlates()
  await coreStore.getSubstrains()
})

const coreStore = useCoreStore()
const selected_plate = ref<string | null>(null)
const selected_substrain = ref<string | null>(null)

const optionsSubstrains = ref<string[]>(coreStore.substrains.map(s => s.name))
const optionsPlates = ref<string[]>(coreStore.plates.map(p => p.barcode))

const {t} = useI18n()

const filterFnPlates = (val: string, update: (fn: () => void) => void) => {
  if (val === '') {
    update(() => {
      optionsPlates.value = coreStore.plates.map(p => p.barcode)
    })
    return
  }
  update(() => {
    const needle = val.toLowerCase()
    optionsPlates.value = coreStore.plates
      .map(p => p.barcode)
      .filter(v => v.toLowerCase().indexOf(needle) > -1)
  })
}

const filterFnSubstrains = (val: string, update: (fn: () => void) => void) => {
  if (val === '') {
    update(() => {
      optionsSubstrains.value = coreStore.substrains.map(s => s.name)
    })
    return
  }
  update(() => {
    const needle = val.toLowerCase()
    optionsSubstrains.value = coreStore.substrains
      .map(s => s.name)
      .filter(v => v.toLowerCase().indexOf(needle) > -1)
  })
}

const onSubmit = async () => {
  if (selected_plate.value) {
    await coreStore.getPlateData(selected_plate.value, selected_substrain.value)
  }
}
</script>

<template>
  <q-card class="my-card">
    <img src="../assets/dna.jpg" alt="virus" style="height: 200px; width: 100%" />

    <q-card-section class="q-pt-none">
      <q-select
        class="tw-my-6"
        color="purple-12"
        v-model="selected_plate"
        use-input
        input-debounce="0"
        label="Plate"
        :options="optionsPlates"
        @filter="filterFnPlates"
        behavior="dialog">
        <template v-slot:no-option>
          <q-item>
            <q-item-section class="text-grey">No results</q-item-section>
          </q-item>
        </template>
      </q-select>
      <q-select
        class="tw-my-6"
        color="purple-12"
        v-model="selected_substrain"
        use-input
        input-debounce="0"
        label="Substrain (optional)"
        :options="optionsSubstrains"
        @filter="filterFnSubstrains"
        behavior="dialog">
        <template v-slot:no-option>
          <q-item>
            <q-item-section class="text-grey">No results</q-item-section>
          </q-item>
        </template>
      </q-select>
    </q-card-section>
    <q-card-actions>
      <q-btn flat color="primary" @click="onSubmit">Submit</q-btn>
    </q-card-actions>
  </q-card>
</template>

<style scoped>
.my-card {
  min-width: 550px;
  margin: auto;
}
</style>
