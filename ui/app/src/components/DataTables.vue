<script setup lang="ts">
import {columnsSampleCount, columnsMetadata, SampleCountsTableColumns} from 'components/data'
import {createSampleCountsData, createMetadataRows} from 'components/helpers'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'
import {computed, ref} from 'vue'
import {useQuasar} from 'quasar'

const {t} = useI18n()
const $q = useQuasar()

const coreStore = useCoreStore()

const selected_sample = ref<string | null>(null)
const optionsSamples = ref<string[]>(coreStore.samples.map(s => s.pseudoanonymized_id))

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
  if (selected_sample.value) {
    coreStore.filterCountDataBySample(selected_sample.value)
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
  <div class="q-pa-md" v-if="coreStore.sampleCounts.length > 0">
    <h3 class="text-h4 text-center q-pb-lg">Counts</h3>

    <div class="tw-mb-10">
      <div class="flex row sample_options_cont">
        <q-select
          class="tw-my-6 sample_select"
          color="purple-12"
          v-model="selected_sample"
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
        <q-btn class="tw-mr-2" color="primary" :label="t('label.apply')" @click="filterBySample"></q-btn>
      </div>
      <q-btn
        icon="query_stats"
        :label="t('label.aggregate')"
        color="primary"
        class="tw-mr-2"
        @click="aggregate"></q-btn>
      <q-btn color="primary" icon="science" :label="t('label.raw_data')" @click="aggregate"></q-btn>
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
    <h3 class="text-h4 text-center q-pb-md">Metadata</h3>
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
.point {
  cursor: pointer;
}
.point:hover {
  background-color: darkslateblue;
  cursor: pointer;
}

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
</style>
