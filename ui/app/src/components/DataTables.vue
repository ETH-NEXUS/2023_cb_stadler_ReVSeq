<script setup lang="ts">
import {columnsSampleCount, columnsMetadata, SampleCountsTableColumns} from 'components/data'
import {createSampleCountsData, createMetadataRows} from 'components/helpers'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'
import {computed, ref} from 'vue'
import {useQuasar} from 'quasar'
import {storeToRefs} from 'pinia'
import PlateFiles from 'components/PlateFiles.vue'

const {t} = useI18n()
const $q = useQuasar()

const coreStore = useCoreStore()

const {selected_sample_id, selected_sample} = storeToRefs(coreStore)
const optionsSamples = ref<string[]>(coreStore.samples.map(s => s.pseudoanonymized_id))
const showPlateFiles = ref<boolean>(false)
const showSampleFiles = ref<boolean>(false)

const aggregate = () => {
  try {
    coreStore.aggregate = !coreStore.aggregate
    coreStore.toggleAggregate()
  } catch (err) {
    console.error(err)
  }
}
const filterColumnNames = computed(() => {
  if (coreStore.aggregate) {
    return columnsSampleCount.filter((c: SampleCountsTableColumns) => c.name !== 'substrain')
  }
  return columnsSampleCount
})
const filter = ref<string>('')

const filterFnSamples = (val: string, update: (fn: () => void) => void) => {
  if (val === '') {
    update(() => {
      optionsSamples.value = coreStore.samples.map(s => s.pseudoanonymized_id)
    })
    return
  }
  update(() => {
    const needle = val.toLowerCase()
    optionsSamples.value = coreStore.samples
      .map(s => s.pseudoanonymized_id)
      .filter(v => v.toLowerCase().indexOf(needle) > -1)
  })
}
const filterBySample = () => {
  if (selected_sample_id.value) {
    coreStore.filterCountDataBySample(selected_sample_id.value)
  } else {
    coreStore.filterCountDataBySample('')
    $q.notify({
      message: 'No sample selected. Resetting the filter.',
      color: 'info',
      icon: 'info',
      position: 'top',
      timeout: 2000,
    })
  }
}
</script>

<template>
  <div class="q-pa-md counts-cont" v-if="coreStore.sampleCounts.length > 0">
    <q-dialog v-model="showPlateFiles">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">{{ t('label.plate_files') }}</div>
          <q-space></q-space>
          <q-btn icon="close" flat round dense v-close-popup></q-btn>
        </q-card-section>

        <q-card-section>
          <PlateFiles domain="plate" />
        </q-card-section>
      </q-card>
    </q-dialog>

    <q-dialog v-model="showSampleFiles">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">{{ t('label.sample_files') }}</div>
          <q-space></q-space>
          <q-btn icon="close" flat round dense v-close-popup></q-btn>
        </q-card-section>

        <q-card-section>
          <PlateFiles domain="sample" />
        </q-card-section>
      </q-card>
    </q-dialog>

    <div class="tw-mb-10">
      <div class="flex row sample_options_cont">
        <q-select
          class="tw-my-6 sample_select"
          color="purple-12"
          v-model="selected_sample_id"
          use-input
          input-debounce="0"
          label="Filter by sample (optional)"
          :options="optionsSamples"
          @filter="filterFnSamples"
          behavior="dialog">
          <template v-slot:no-option>
            <q-item>
              <q-item-section class="text-grey">No results</q-item-section>
            </q-item>
          </template>
        </q-select>
        <q-btn
          class="tw-mr-2 tw-mt-4"
          color="primary"
          :label="t('label.apply')"
          @click="filterBySample"></q-btn>
      </div>
      <div class="flex row wrap">
        <q-btn
          icon="query_stats"
          :label="t('label.aggregate')"
          color="primary"
          class="tw-mr-2 tw-my-2"
          @click="aggregate"></q-btn>
        <q-btn
          class="tw-my-2 tw-mr-2"
          color="primary"
          icon="science"
          :label="t('label.raw_data')"
          @click="aggregate"></q-btn>
        <q-btn
          icon="file_download"
          class="tw-my-2 tw-mr-2"
          :label="t('label.plate_files')"
          color="primary"
          @click="showPlateFiles = true"></q-btn>

        <q-btn
          v-if="selected_sample"
          icon="file_download"
          class="tw-my-2"
          :label="t('label.sample_files')"
          color="primary"
          @click="showSampleFiles = true"></q-btn>
      </div>
    </div>
    <q-table
      :filter="filter"
      :rows-per-page-options="[10, 20, 50, 100]"
      rows-per-page="20"
      flat
      bordered
      :rows="createSampleCountsData(coreStore.tableData)"
      :columns="filterColumnNames"
      row-key="name"></q-table>
  </div>
  <div class="q-pa-md" v-if="coreStore.metadata.length > 0">
    <h3 class="text-h4 text-center q-pb-md">{{ t('titles.metadata') }}</h3>
    <q-table
      :rows-per-page-options="[50, 100]"
      rows-per-page="50"
      flat
      bordered
      :rows="createMetadataRows(coreStore.metadata)"
      :columns="columnsMetadata"
      row-key="name"></q-table>
  </div>
</template>

<style scoped>
.sample_options_cont {
  width: 40%;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
}

.sample_select {
  flex: 1;
  margin-right: 14px;
}

q-btn {
  flex-shrink: 0; /* Prevents the button from shrinking if space is tight */
}

.counts-cont {
  max-width: 80%;
  margin-bottom: 5px;
}
</style>
